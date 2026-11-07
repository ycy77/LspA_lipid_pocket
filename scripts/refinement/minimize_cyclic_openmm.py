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
    # Environment switches: hydrogen addition defaults off; first-residue terminal-H cleanup defaults on.
    add_hydrogens = os.environ.get('ADD_HYDROGENS', '0').lower() in ['1','true','yes']
    strip_terminal_h = os.environ.get('STRIP_TERMINAL_H', '1').lower() in ['1','true','yes']
    try:
        # 1. Repair missing atoms with PDBFixer.
        fixer = PDBFixer(filename=input_pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()  # Add missing heavy atoms; hydrogen addition is handled separately.
        fixer.removeHeterogens(False)  # Remove heterogens and retain the protein.
        modeller = Modeller(fixer.topology, fixer.positions)
        target_chain = None

        # 2. Locate the target chain.
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

                # 3. Remove OXT.
                atoms_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                if atoms_to_delete:
                    modeller.delete(atoms_to_delete)
                    print(f"Deleted OXT from {last_res.name} {last_res.id}")

                # 4. Add the ring-closing peptide bond if absent.
                n_atom = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                c_atom = next((a for a in last_res.atoms() if a.name.strip().upper() == 'C'), None)

                if n_atom and c_atom:
                    # Print the N-C distance for inspection.
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

                    # Check whether the peptide bond already exists.
                    existing_bonds = list(modeller.topology.bonds())
                    bond_exists = any((b.atom1 == n_atom and b.atom2 == c_atom) or
                                      (b.atom1 == c_atom and b.atom2 == n_atom) for b in existing_bonds)
                    if not bond_exists:
                        modeller.topology.addBond(n_atom, c_atom)
                        print(f"Added peptide bond between {first_res.name} N and {last_res.name} C")
                    else:
                        print(f"Peptide bond already exists for chain {cyclize_chain_id}, skipping addBond")

                    # Optionally remove first-residue terminal hydrogens to avoid N-terminal Amber template assignment.
                    if strip_terminal_h:
                        h_names = {"H", "H1", "H2", "H3", "HN", "HT1", "HT2", "HT3"}
                        h_to_delete = [a for a in first_res.atoms() if a.name.strip().upper() in h_names]
                        if h_to_delete:
                            modeller.delete(h_to_delete)
                            print(f"Stripped terminal H from first residue {first_res.name} {first_res.id}")
                else:
                    print(f"Warning: N/C atoms not found in chain {cyclize_chain_id}")

        # 5. Load the protein force field.
        forcefield = ForceField(
            'amber14/protein.ff14SB.xml',
            'amber14/tip3p.xml'
        )

        # Make hydrogen addition optional (off by default) to avoid rebuilding incompatible terminal groups.
        if add_hydrogens:
            try:
                modeller.addHydrogens(forcefield, pH=7.0)
                print("Hydrogens added using Modeller (amber14)")
            except Exception as e:
                print(f"Warning: addHydrogens failed: {e}")
        else:
            print("Hydrogen addition skipped by env ADD_HYDROGENS=0")

        # 6. Create the system, with a hydrogen-addition fallback for missing templates.
        # Choose the nonbonded method based on periodic box vectors in the topology.
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
                    # Second fallback: add hydrogens to the current structure using PDBFixer.
                    try:
                        buf = io.StringIO()
                        # Preserve chain and residue identifiers for locating the target chain later.
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

                # Normalize termini after cyclization: remove OXT and normalize N-terminal hydrogens.
                # After the PDBFixer fallback, locate the target chain again.
                # Retrieve the target chain again.
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
                        # Retrieve residues from the updated chain reference.
                        if used_fix_method == "pdbfixer":
                            residues = list(target_chain_ref.residues())
                        first_res = residues[0]
                        last_res = residues[-1]
                        # Restore the ring-closing bond if PDBFixer removed the inter-residue bond.
                        n_atom_fb = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                        c_atom_fb = next((a for a in last_res.atoms() if a.name.strip().upper() == 'C'), None)
                        if n_atom_fb and c_atom_fb:
                            existing_bonds_fb = list(modeller.topology.bonds())
                            bond_exists_fb = any((b.atom1 == n_atom_fb and b.atom2 == c_atom_fb) or
                                                 (b.atom1 == c_atom_fb and b.atom2 == n_atom_fb) for b in existing_bonds_fb)
                            if not bond_exists_fb:
                                modeller.topology.addBond(n_atom_fb, c_atom_fb)
                                print(f"Fallback normalization: re-added peptide bond {first_res.name} N — {last_res.name} C")
                        # Remove terminal OXT if atom or hydrogen addition restored it.
                        oxt_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                        if oxt_to_delete:
                            modeller.delete(oxt_to_delete)
                            print(f"Fallback normalization: removed OXT from {last_res.name} {last_res.id}")

                        # Normalize first-residue N hydrogens: none for PRO, at most one for other residues.
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

                # Retry system creation after the fallback, checking periodic box vectors again.
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
                # Update the target chain reference after successful system creation for subsequent restraints.
                if target_chain_ref is not None:
                    target_chain = target_chain_ref
            else:
                # Reraise errors unrelated to hydrogen assignment.
                raise

        # Normalize termini after cyclization: remove OXT and normalize N-terminal hydrogens.
        if target_chain is not None:
            residues = list(target_chain.residues())
            if len(residues) >= 2:
                first_res = residues[0]
                last_res = residues[-1]
                # Remove terminal OXT if atom or hydrogen addition restored it.
                oxt_to_delete = [a for a in last_res.atoms() if a.name.strip().upper() == 'OXT']
                if oxt_to_delete:
                    modeller.delete(oxt_to_delete)
                    print(f"Terminal normalization: removed OXT from {last_res.name} {last_res.id} (post-fix/hydrogens)")

                # Normalize first-residue N hydrogens: none for PRO, at most one for other residues.
                n_atom = next((a for a in first_res.atoms() if a.name.strip().upper() == 'N'), None)
                if n_atom is not None:
                    # Find hydrogens bonded to N.
                    bonded = []
                    for b in modeller.topology.bonds():
                        if b.atom1 is n_atom:
                            bonded.append(b.atom2)
                        elif b.atom2 is n_atom:
                            bonded.append(b.atom1)
                    hydrogens_on_n = [a for a in bonded if getattr(a, 'element', None) is not None and a.element is app.Element.getByAtomicNumber(1)]

                    if first_res.name.strip().upper().startswith('PRO'):
                        # PRO nitrogen should have no attached hydrogen in this cyclic topology.
                        if hydrogens_on_n:
                            modeller.delete(hydrogens_on_n)
                            print(f"Terminal normalization: removed all N-H on PRO {first_res.id}")
                    else:
                        # Keep at most one N-H for other residues, preferring H or HN.
                        keep = None
                        for cand_name in ('H', 'HN'):
                            keep = next((h for h in hydrogens_on_n if h.name.strip().upper() == cand_name), None)
                            if keep:
                                break
                        # If no preferred name is found, retain the first hydrogen.
                        if hydrogens_on_n:
                            to_delete = [h for h in hydrogens_on_n if h is not keep]
                            if len(hydrogens_on_n) > 1 and to_delete:
                                modeller.delete(to_delete)
                                print(f"Terminal normalization: pruned extra N-H on {first_res.name} {first_res.id}, kept {keep.name if keep else 'one'}")

        # Optional N-O distance restraint to reduce overlap at the cyclization site.
        # Parameters: enable switch, target r0 (nm), spring constant k (kJ/mol/nm^2), and one-sided threshold.
        # Enabled by default; supports a one-sided flat-bottom restraint active only below r_min.
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
                                    # One-sided flat bottom: repel only when r < r0, with spring constant k_no.
                                    expr = 'step(r0-r)*k_no*(r0-r)^2'
                                else:
                                    # Two-sided harmonic restraint.
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

        # 7. Fix chains other than the cyclic target chain.
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

        # 8. Optional local relaxation: release terminal residues and their neighbors while restraining the rest.
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

        # Keep the standard minimization procedure without the former O-kick coordinate adjustment.

        # 9. Set up simulation and minimize energy.
        integrator = openmm.LangevinIntegrator(
            300*unit.kelvin,
            1.0/unit.picoseconds,
            0.002*unit.picoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        # Support CUDA and CPU platforms.
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
                # Fall back to the default platform.
                simulation = Simulation(modeller.topology, system, integrator)
                print("CUDA/CPU not available, using default platform")

        simulation.context.setPositions(modeller.positions)
        # Validate the initial energy.
        state = simulation.context.getState(getEnergy=True)
        try:
            print(f"Initial energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
        except Exception:
            print("Initial energy: (unavailable)")
        # Minimize in two stages: local relaxation, then global minimization.
        if warmup_force is not None:
            try:
                print(f"Warmup minimization: {warmup_iters} iterations with local restraints")
                simulation.minimizeEnergy(maxIterations=warmup_iters)
            except Exception as e:
                print(f"Warmup minimization failed: {e}")
            # Disable local restraints by setting their global parameter to zero.
            try:
                simulation.context.setParameter('k_warmup', 0.0)
                print("Warmup restraints disabled (k_warmup=0)")
            except Exception:
                print("Failed to disable warmup via global parameter; continuing")

        # Stage 2: global minimization.
        simulation.minimizeEnergy(maxIterations=10000)
        # Record the final energy after stage 2.
        try:
            state_final = simulation.context.getState(getEnergy=True)
            print(f"Final energy: {state_final.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.2f} kJ/mol")
        except Exception:
            print("Final energy: (unavailable)")

        # 10. Save the result.
        # Use minimized coordinates and remove all hydrogens before saving.
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
