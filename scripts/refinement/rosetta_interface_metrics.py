
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path

import os
import csv
import multiprocessing
from multiprocessing import Pool
import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap

# 每个进程的初始化函数
def init_worker():
    # 每个进程需要单独初始化Pyrosetta
    pyrosetta.init()

# 计算单个结构的指标
def process_pdb(args):
    backbone, best_model, pdb_path = args
    try:
        # 加载结构
        pose = pyrosetta.Pose()
        pyrosetta.pose_from_file(pose, pdb_path)

        # 创建得分函数
        scorefxn = ScoreFunctionFactory.create_score_function("ref2015")

        # 创建MoveMap，允许所有原子移动
        movemap = MoveMap()
        movemap.set_bb(True)  # 允许 backbone 移动
        movemap.set_chi(True)  # 允许侧链移动
        movemap.set_jump(True)  # 允许链间移动

        # 能量最小化
        min_mover = MinMover()
        min_mover.movemap(movemap)
        min_mover.score_function(scorefxn)
        min_mover.min_type("dfpmin_armijo_nonmonotone")
        min_mover.tolerance(0.1)
        min_mover.max_iter(1000)
        min_mover.apply(pose)

        # 创建InterfaceAnalyzerMover来计算界面指标
        interface_mover = InterfaceAnalyzerMover()
        interface_mover.set_scorefunction(scorefxn)
        interface_mover.set_interface("A_B")  # 假设界面是A链和B链之间
        interface_mover.apply(pose)

        # 获取dSASA (界面溶剂可及表面积变化)
        dsasa = interface_mover.get_interface_delta_sasa()

        # 获取ddG (界面能量变化)
        ddg = interface_mover.get_interface_dG()

        # 获取BUH (埋藏未满足氢键数)
        buh = interface_mover.get_interface_delta_hbond_unsat()

        # Legacy CSV name CYCttE: this is the total ref2015 score, not a cycle-specific energy.
        cyctte = scorefxn(pose)

        return (backbone, best_model, dsasa, ddg, buh, cyctte, True, "")
    except Exception as e:
        error_msg = str(e)
        return (backbone, best_model, "N/A", "N/A", "N/A", "N/A", False, error_msg)

# 主函数
def main():
    # 读取final_summary_with_Bcount.tsv文件
    tsv_file = env_path("LSP_SUMMARY_TSV", f"{DATA}/cyclic_youhua_1/final_summary.tsv")
    pdb_files = []

    with open(tsv_file, "r", encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            backbone = row['Backbone']
            best_model = row['Best_Model']
            # 构建完整路径
            pdb_path = os.path.join(env_path("LSP_PDB_ROOT", f"{DATA}/cyclic_youhua_1"), backbone, os.environ.get("LSP_MODEL_SUBDIR", "minimized"), best_model)
            if os.path.exists(pdb_path):
                pdb_files.append((backbone, best_model, pdb_path))

    print(f"找到 {len(pdb_files)} 个PDB文件")

    # 定义输出文件
    output_file = env_path("LSP_OUTPUT_CSV", f"{ANALYSIS}/pyrosetta_youhua_results_minimized2.csv")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if not pdb_files:
        raise SystemExit("No structures matched the summary; check LSP_PDB_ROOT and LSP_MODEL_SUBDIR")

    # 写入表头
    with open(output_file, "w") as f:
        f.write("Backbone,Best_Model,dSASA,ddG,BUH,CYCttE\n")

    # 准备参数
    params = [(backbone, best_model, pdb_path) for backbone, best_model, pdb_path in pdb_files]

    # 使用多进程处理
    num_processes = max(1, min(multiprocessing.cpu_count(), int(os.environ.get("LSP_WORKERS", "4"))))  # 使用最多36个进程
    print(f"使用 {num_processes} 个进程进行计算")

    with Pool(processes=num_processes, initializer=init_worker) as pool:
        results = pool.imap(process_pdb, params, chunksize=1)

        # 处理结果
        for i, result in enumerate(results):
            backbone, best_model, dsasa, ddg, buh, cyctte, success, error_msg = result
            progress = (i + 1) / len(pdb_files) * 100

            if success:
                print(f"[{progress:.1f}%] 成功: {backbone} - dSASA={dsasa:.2f}, ddG={ddg:.2f}, BUH={buh}, CYCttE={cyctte:.2f}")
                # 写入成功结果
                with open(output_file, "a") as f:
                    f.write(f"{backbone},{best_model},{dsasa},{ddg},{buh},{cyctte}\n")
            else:
                print(f"[{progress:.1f}%] 失败: {backbone} - {error_msg}")
                # 写入失败结果
                with open(output_file, "a") as f:
                    f.write(f"{backbone},{best_model},N/A,N/A,N/A,N/A\n")

    print(f"\n分析完成，结果已保存到 {output_file}")
    return output_file

if __name__ == "__main__":
    main()
