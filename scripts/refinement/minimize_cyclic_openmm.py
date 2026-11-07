#!/usr/bin/env python3

# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path

import os
import io
import sys
import multiprocessing
from openmm import app, openmm, unit
from openmm.app import Simulation, PDBFile, ForceField
from pdbfixer import PDBFixer
from openmm.app import Modeller

def minimize_cyclic_pdb(input_pdb_path, output_pdb_path, cyclize_chain_id='B'):
    print(f"Processing {input_pdb_path} -> {output_pdb_path}")
    # 环境开关：是否补氢，默认跳过；是否清理首残终端氢，默认开启
    add_hydrogens = os.environ.get('ADD_HYDROGENS', '0').lower() in ['1','true','yes']
    strip_terminal_h = os.environ.get('STRIP_TERMINAL_H', '1').lower() in ['1','true','yes']
    try:
        # 1. 用PDBFixer彻底修复PDB（补全原子）
        fixer = PDBFixer(filename=input_pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()  # 补全缺失的氢和重原子
        fixer.removeHeterogens(False)  # 保留蛋白质
        modeller = Modeller(fixer.topology, fixer.positions)
        target_chain = None

        # 2. 查找目标链
        for ch in modeller.topology.chains():
            if str(ch.id).strip() == cyclize_chain_id.strip():
                target_chain = ch
                break

        if target_chain is None:
            print(f"Warning: chain '{cyclize_chain_id}' not found. Skipping cyclization.")
        else:
            residues = list(target_chain.residues())
            if len(residues) >= 2:
                first_res = residues[0]
                last_res = residues[-1]
                print(f"Cyclizing chain B: first={first_res.name} {first_res.id}, last={last_res.name} {last_res.id}")

                # 3. 删除OXT
                atoms_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                if atoms_to_delete:
                    modeller.delete(atoms_to_delete)
                    print(f"Deleted OXT from {last_res.name} {last_res.id}")

                # 4. 添加闭环肽键（避免重复加键）
                n_atom = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                c_atom = next((a for a in last_res.atoms() if a.name.strip().upper() == 'C'), None)

                if n_atom and c_atom:
                    # 打印 N-C 距离用于验证
                    atom_list = list(modeller.topology.atoms())
                    atom_index_map = {a: i for i, a in enumerate(atom_list)}
                    n_index = atom_index_map.get(n_atom)
                    c_index = atom_index_map.get(c_atom)
                    if n_index is not None and c_index is not None:
                        try:
                            nc_dist_nm = unit.norm(modeller.positions[n_index] - modeller.positions[c_index]).value_in_unit(unit.nanometers)
                            print(f"Topology: N(first)-C(last) geometric distance ≈ {nc_dist_nm:.3f} nm for chain {cyclize_chain_id}.")
                        except Exception:
                            pass

                    # 检测是否已有肽键
                    existing_bonds = list(modeller.topology.bonds())
                    bond_exists = any((b.atom1 == n_atom and b.atom2 == c_atom) or
                                      (b.atom1 == c_atom and b.atom2 == n_atom) for b in existing_bonds)
                    if not bond_exists:
                        modeller.topology.addBond(n_atom, c_atom)
                        print(f"Added peptide bond between {first_res.name} N and {last_res.name} C")
                    else:
                        print(f"Peptide bond already exists for chain {cyclize_chain_id}, skipping addBond")

                    # 可选：清理首残终端氢，避免 Amber 模板误识别为 N 端
                    if strip_terminal_h:
                        h_names = {"H", "H1", "H2", "H3", "HN", "HT1", "HT2", "HT3"}
                        h_to_delete = [a for a in first_res.atoms() if a.name.strip().upper() in h_names]
                        if h_to_delete:
                            modeller.delete(h_to_delete)
                            print(f"Stripped terminal H from first residue {first_res.name} {first_res.id}")
                else:
                    print(f"Warning: N/C atoms not found in chain {cyclize_chain_id}")

        # 5. 加载蛋白质专用力场
        forcefield = ForceField(
            'amber14/protein.ff14SB.xml',
            'amber14/tip3p.xml'
        )

        # 关键修正：补氢可控（默认跳过），避免被补氢重建端基导致模板不匹配
        if add_hydrogens:
            try:
                modeller.addHydrogens(forcefield, pH=7.0)
                print("Hydrogens added using Modeller (amber14)")
            except Exception as e:
                print(f"Warning: addHydrogens failed: {e}")
        else:
            print("Hydrogen addition skipped by env ADD_HYDROGENS=0")

        # 6. 创建系统（带模板缺失的自动补氢回退）
        # 根据拓扑是否包含周期盒向量选择非键相互作用方法
        has_box = None
        try:
            has_box = modeller.topology.getPeriodicBoxVectors()
        except Exception:
            has_box = None
        if has_box is None:
            nb_method = app.CutoffNonPeriodic
            print("Detected non-periodic topology; using CutoffNonPeriodic for nonbonded interactions")
        else:
            nb_method = app.PME
            print("Periodic box detected; using PME for nonbonded interactions")
        try:
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=nb_method,
                nonbondedCutoff=1.0*unit.nanometers,
                constraints=app.HBonds,
                rigidWater=True,
                ignoreExternalBonds=True
            )
        except Exception as e:
            msg = str(e)
            print(f"ForceField system build failed: {msg}")
            lower = msg.lower()
            needs_h_fallback = ("template" in lower) or ("missing" in lower and "hydrogen" in lower)
            if needs_h_fallback:
                print("Triggering fallback: adding hydrogens and re-normalizing terminals")
                used_fix_method = None
                try:
                    modeller.addHydrogens(forcefield, pH=7.0)
                    print("Fallback: hydrogens added (amber14)")
                    used_fix_method = "modeller"
                except Exception as e2:
                    print(f"Fallback addHydrogens failed: {e2}")
                    # 二级回退：用 PDBFixer 在当前结构上补氢
                    try:
                        buf = io.StringIO()
                        # 保留原始链/残基编号以便后续定位目标链
                        PDBFile.writeFile(modeller.topology, modeller.positions, buf, keepIds=True)
                        pdb_str = buf.getvalue()
                        fixer2 = PDBFixer(pdbfile=io.StringIO(pdb_str))
                        fixer2.findMissingResidues()
                        fixer2.findMissingAtoms()
                        fixer2.addMissingAtoms()
                        fixer2.addMissingHydrogens()
                        modeller = Modeller(fixer2.topology, fixer2.positions)
                        print("Fallback: hydrogens added using PDBFixer")
                        used_fix_method = "pdbfixer"
                    except Exception as e3:
                        print(f"Secondary fallback (PDBFixer) failed: {e3}")

                # 再次规范端基：确保闭环后按内部残基处理（删除 OXT、规整 N 端氢）
                # 注意：如使用PDBFixer回退，需重新定位目标链对象
                # 先重新获取目标链
                target_chain_ref = None
                chain_ids = []
                for ch in modeller.topology.chains():
                    if str(ch.id).strip() == cyclize_chain_id.strip():
                        target_chain_ref = ch
                        break
                    chain_ids.append(str(ch.id).strip())
                if target_chain_ref is None:
                    print(f"Fallback normalization: target chain '{cyclize_chain_id}' not found; available chains: {chain_ids}")

                if target_chain_ref is not None:
                    residues = list(target_chain.residues())
                    if len(residues) >= 2:
                        # 若使用新的链引用，重新取残基
                        if used_fix_method == "pdbfixer":
                            residues = list(target_chain_ref.residues())
                        first_res = residues[0]
                        last_res = residues[-1]
                        # 重新确保闭环肽键存在（PDBFixer可能丢失跨残基键）
                        n_atom_fb = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                        c_atom_fb = next((a for a in last_res.atoms() if a.name.strip().upper() == 'C'), None)
                        if n_atom_fb and c_atom_fb:
                            existing_bonds_fb = list(modeller.topology.bonds())
                            bond_exists_fb = any((b.atom1 == n_atom_fb and b.atom2 == c_atom_fb) or
                                                 (b.atom1 == c_atom_fb and b.atom2 == n_atom_fb) for b in existing_bonds_fb)
                            if not bond_exists_fb:
                                modeller.topology.addBond(n_atom_fb, c_atom_fb)
                                print(f"Fallback normalization: re-added peptide bond {first_res.name} N — {last_res.name} C")
                        # 删除末端 OXT（若在补氢或补原子后被重新添加）
                        oxt_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                        if oxt_to_delete:
                            modeller.delete(oxt_to_delete)
                            print(f"Fallback normalization: removed OXT from {last_res.name} {last_res.id}")

                        # 规整首残 N 上的氢：PRO 不应有氢，其它残基最多保留一个
                        n_atom = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                        if n_atom is not None:
                            bonded = []
                            for b in modeller.topology.bonds():
                                if b.atom1 is n_atom:
                                    bonded.append(b.atom2)
                                elif b.atom2 is n_atom:
                                    bonded.append(b.atom1)
                            hydrogens_on_n = [a for a in bonded if getattr(a, 'element', None) is not None and a.element is app.Element.getByAtomicNumber(1)]

                            if first_res.name.strip().upper().startswith('PRO'):
                                if hydrogens_on_n:
                                    modeller.delete(hydrogens_on_n)
                                    print(f"Fallback normalization: removed all N-H on PRO {first_res.id}")
                            else:
                                keep = None
                                for cand_name in ('H', 'HN'):
                                    keep = next((h for h in hydrogens_on_n if h.name.strip().upper() == cand_name), None)
                                    if keep:
                                        break
                                if hydrogens_on_n:
                                    to_delete = [h for h in hydrogens_on_n if h is not keep]
                                    if len(hydrogens_on_n) > 1 and to_delete:
                                        modeller.delete(to_delete)
                                        print(f"Fallback normalization: pruned extra N-H on {first_res.name} {first_res.id}, kept {keep.name if keep else 'one'}")

                # 回退后重试创建系统（重新检测周期盒）
                has_box_fb = None
                try:
                    has_box_fb = modeller.topology.getPeriodicBoxVectors()
                except Exception:
                    has_box_fb = None
                if has_box_fb is None:
                    nb_method_fb = app.CutoffNonPeriodic
                    print("Fallback: non-periodic topology; using CutoffNonPeriodic")
                else:
                    nb_method_fb = app.PME
                    print("Fallback: periodic box detected; using PME")
                try:
                    system = forcefield.createSystem(
                        modeller.topology,
                        nonbondedMethod=nb_method_fb,
                        nonbondedCutoff=1.0*unit.nanometers,
                        constraints=app.HBonds,
                        rigidWater=True,
                        ignoreExternalBonds=True
                    )
                    print("Fallback: system build succeeded after adding hydrogens")
                except Exception as e3:
                    print(f"Fallback: system build still failing: {e3}")
                    raise
                # 如果成功构建系统，更新目标链引用供后续约束使用
                if target_chain_ref is not None:
                    target_chain = target_chain_ref
            else:
                # 与氢无关的错误，直接抛出
                raise

        # 再次规范端基：确保闭环后按内部残基处理（删除 OXT、规整 N 端氢）
        if target_chain is not None:
            residues = list(target_chain.residues())
            if len(residues) >= 2:
                first_res = residues[0]
                last_res = residues[-1]
                # 删除末端 OXT（若在补氢或补原子后被重新添加）
                oxt_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                if oxt_to_delete:
                    modeller.delete(oxt_to_delete)
                    print(f"Terminal normalization: removed OXT from {last_res.name} {last_res.id} (post-fix/hydrogens)")

                # 规整首残 N 上的氢：PRO 不应有氢，其它残基最多保留一个
                n_atom = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                if n_atom is not None:
                    # 找到与 N 相连的氢
                    bonded = []
                    for b in modeller.topology.bonds():
                        if b.atom1 is n_atom:
                            bonded.append(b.atom2)
                        elif b.atom2 is n_atom:
                            bonded.append(b.atom1)
                    hydrogens_on_n = [a for a in bonded if getattr(a, 'element', None) is not None and a.element is app.Element.getByAtomicNumber(1)]

                    if first_res.name.strip().upper().startswith('PRO'):
                        # PRO 的 N 不应带氢
                        if hydrogens_on_n:
                            modeller.delete(hydrogens_on_n)
                            print(f"Terminal normalization: removed all N-H on PRO {first_res.id}")
                    else:
                        # 其它残基最多保留一个 N-H（优先保留名字为 H/HN 的）
                        keep = None
                        for cand_name in ('H', 'HN'):
                            keep = next((h for h in hydrogens_on_n if h.name.strip().upper() == cand_name), None)
                            if keep:
                                break
                        # 若未找到命名合适者，则保留第一个
                        if hydrogens_on_n:
                            to_delete = [h for h in hydrogens_on_n if h is not keep]
                            if len(hydrogens_on_n) > 1 and to_delete:
                                modeller.delete(to_delete)
                                print(f"Terminal normalization: pruned extra N-H on {first_res.name} {first_res.id}, kept {keep.name if keep else 'one'}")

        # 可选：N–O 距离约束以缓解环化位点重叠
        # 参数：开启开关、目标距离 r0（nm）、弹簧常数 k（kJ/mol/nm^2）、是否仅在“过近”时施加以及判定阈值
        # 默认开启；支持平底单侧模式（只在 r < r_min 时施加）
        no_restraint = os.environ.get('NO_DIST_RESTRAINT', '1').lower() in ['1','true','yes']
        no_r0_nm = float(os.environ.get('NO_R0_NM', '0.30'))
        no_k = float(os.environ.get('NO_K', '6000.0'))
        no_only_when_close = os.environ.get('NO_ONLY_WHEN_CLOSE', '0').lower() in ['1','true','yes']
        no_close_th_nm = float(os.environ.get('NO_CLOSE_TH_NM', '0.26'))
        no_flat_bottom = os.environ.get('NO_FLAT_BOTTOM', '1').lower() in ['1','true','yes']
        if target_chain is not None and no_restraint:
            residues = list(target_chain.residues())
            if len(residues) >= 2:
                first_res = residues[0]
                last_res = residues[-1]
                n_atom_no = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                o_atom_no = next((a for a in last_res.atoms() if a.name.strip().upper() == 'O'), None)
                if n_atom_no and o_atom_no:
                    try:
                        atom_list = list(modeller.topology.atoms())
                        atom_index_map = {a: i for i, a in enumerate(atom_list)}
                        n_idx = atom_index_map.get(n_atom_no)
                        o_idx = atom_index_map.get(o_atom_no)
                        if n_idx is not None and o_idx is not None:
                            r_no = unit.norm(modeller.positions[n_idx] - modeller.positions[o_idx]).value_in_unit(unit.nanometers)
                            apply_this = True
                            if no_only_when_close and not (r_no < no_close_th_nm):
                                apply_this = False
                            if apply_this:
                                if no_flat_bottom:
                                    # 单侧平底：只在 r < r0 时推开，力常数 k_no
                                    expr = 'step(r0-r)*k_no*(r0-r)^2'
                                else:
                                    # 双侧谐波
                                    expr = 'k_no*(r-r0)^2'
                                no_force = openmm.CustomBondForce(expr)
                                no_force.addGlobalParameter('k_no', no_k)
                                no_force.addGlobalParameter('r0', no_r0_nm)
                                no_force.addBond(n_idx, o_idx, [])
                                system.addForce(no_force)
                                mode = 'flat-bottom (one-sided)' if no_flat_bottom else 'harmonic'
                                print(f"N–O restraint added ({mode}): r_init={r_no:.3f} nm, r0={no_r0_nm:.3f} nm, k={no_k}")
                            else:
                                print(f"N–O initial distance {r_no:.3f} nm not close; no restraint applied")
                    except Exception as e:
                        print(f"Failed to add N–O restraint: {e}")

        # 7. 固定非闭环链
        if target_chain is not None:
            force = openmm.CustomExternalForce('1000*(x-x0)^2 + 1000*(y-y0)^2 + 1000*(z-z0)^2')
            force.addPerParticleParameter('x0')
            force.addPerParticleParameter('y0')
            force.addPerParticleParameter('z0')
            for i, atom in enumerate(modeller.topology.atoms()):
                if atom.residue.chain.id != target_chain.id:
                    force.addParticle(i, modeller.positions[i])
            system.addForce(force)
            print(f"Fixed all chains except '{cyclize_chain_id}'")
        else:
            print("No target chain found, skipping position restraints")

        # 8.（可选）闭环局部预约束：仅释放首/末残及其邻居，其余目标链原位约束，先做局部松弛
        local_warmup = os.environ.get('LOCAL_WARMUP', '1').lower() in ['1','true','yes']
        warmup_iters = int(os.environ.get('LOCAL_WARMUP_ITERS', '3000'))
        warmup_k = float(os.environ.get('LOCAL_WARMUP_K', '2000.0'))  # kJ/mol/nm^2
        warmup_force = None
        if target_chain is not None and local_warmup:
            residues_all = list(target_chain.residues())
            exclude_idx = {0, len(residues_all)-1}
            if len(residues_all) >= 3:
                exclude_idx.update({1, len(residues_all)-2})
            exclude_set = {residues_all[i] for i in exclude_idx}

            warmup_force = openmm.CustomExternalForce('k_warmup*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
            warmup_force.addPerParticleParameter('x0')
            warmup_force.addPerParticleParameter('y0')
            warmup_force.addPerParticleParameter('z0')
            warmup_force.addGlobalParameter('k_warmup', warmup_k)
            added_count = 0
            for i, atom in enumerate(modeller.topology.atoms()):
                if atom.residue.chain.id == target_chain.id and atom.residue not in exclude_set:
                    warmup_force.addParticle(i, modeller.positions[i])
                    added_count += 1
            if added_count > 0:
                system.addForce(warmup_force)
                print(f"Local warmup restraints added on target chain (atoms: {added_count}, k={warmup_k}): free residues idx {sorted(exclude_idx)}")
            else:
                print("Local warmup requested, but no atoms qualified; skipping warmup stage")

        #（回退）移除 O-kick 坐标微调，保持原始最小化流程

        # 9. 模拟与能量最小化
        integrator = openmm.LangevinIntegrator(
            300*unit.kelvin,
            1.0/unit.picoseconds,
            0.002*unit.picoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        # 兼容CUDA/CPU
        try:
            platform = openmm.Platform.getPlatformByName('CUDA')
            properties = {'Precision': 'mixed'}
            simulation = Simulation(modeller.topology, system, integrator, platform, properties)
        except Exception:
            try:
                platform = openmm.Platform.getPlatformByName('CPU')
                simulation = Simulation(modeller.topology, system, integrator, platform)
                print("CUDA not available, using CPU platform")
            except Exception:
                # 最后回退到默认平台
                simulation = Simulation(modeller.topology, system, integrator)
                print("CUDA/CPU not available, using default platform")

        simulation.context.setPositions(modeller.positions)
        # 验证初始能量
        state = simulation.context.getState(getEnergy=True)
        try:
            print(f"Initial energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
        except Exception:
            print("Initial energy: (unavailable)")
        # 能量最小化（两阶段：先局部，再全局）
        if warmup_force is not None:
            try:
                print(f"Warmup minimization: {warmup_iters} iterations with local restraints")
                simulation.minimizeEnergy(maxIterations=warmup_iters)
            except Exception as e:
                print(f"Warmup minimization failed: {e}")
            # 关闭局部约束（将全局参数置零）
            try:
                simulation.context.setParameter('k_warmup', 0.0)
                print("Warmup restraints disabled (k_warmup=0)")
            except Exception:
                print("Failed to disable warmup via global parameter; continuing")

        # 阶段 2：全局最小化
        simulation.minimizeEnergy(maxIterations=10000)
        # 记录阶段 2 最终能量
        try:
            state_final = simulation.context.getState(getEnergy=True)
            print(f"Final energy: {state_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
        except Exception:
            print("Final energy: (unavailable)")

        # 10. 保存结果
        # 使用最小化后的坐标，并在保存前移除所有氢原子
        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions()
        modeller_out = Modeller(modeller.topology, positions)
        hydrogen = app.Element.getByAtomicNumber(1)
        h_to_delete = [a for a in modeller_out.topology.atoms() if getattr(a, 'element', None) is hydrogen]
        if h_to_delete:
            modeller_out.delete(h_to_delete)
            print(f"Stripped all hydrogens before saving (removed {len(h_to_delete)} atoms)")

        os.makedirs(os.path.dirname(output_pdb_path), exist_ok=True)
        with open(output_pdb_path, 'w') as f:
            PDBFile.writeFile(modeller_out.topology, modeller_out.positions, f, keepIds=True)
        print(f"Successfully saved to {output_pdb_path}\n")

    except Exception as e:
        print(f"Error processing {input_pdb_path}: {str(e)}", file=sys.stderr)

# Run sequentially: avoids Queue.empty() races and makes failures visible.
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Original OpenMM cyclic-peptide minimization")
    parser.add_argument("--input-dir", required=True, type=_Path)
    parser.add_argument("--output-dir", required=True, type=_Path)
    parser.add_argument("--chain", default="B")
    args = parser.parse_args()
    inp, out = args.input_dir.resolve(), args.output_dir.resolve()
    if out == inp or inp in out.parents:
        parser.error("Output must be outside the input tree")
    pdbs = sorted(inp.rglob("*.pdb"))
    if not pdbs:
        parser.error(f"No PDB files under {inp}")
    failures = 0
    for pdb in pdbs:
        dest = out / pdb.relative_to(inp).parent / "minimized" / (pdb.stem + "_minimized.pdb")
        if dest.is_file() and dest.stat().st_size:
            print(f"SKIP {dest}")
            continue
        minimize_cyclic_pdb(str(pdb), str(dest), args.chain)
        if not dest.is_file() or not dest.stat().st_size:
            failures += 1
    print(f"Processed {len(pdbs)} structures; failures={failures}")
    raise SystemExit(1 if failures else 0)
