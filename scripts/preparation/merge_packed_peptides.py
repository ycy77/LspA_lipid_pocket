#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：
批量处理侧链包装PDB与原始B链的对齐替换：
1. 将每个packed PDB的B链与原始PDB的B链对齐
2. 替换原始PDB的B链为对齐后的packed B链
3. 结果保存至每个样本的merge_pdb目录
"""


# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path

import os
import glob
import pymol
from pymol import cmd
import logging

# ============ 1. 路径配置 ============
BASE_DIR = f"{WORKFLOW}"
# 原始PDB所在目录（含待对齐的B链）
ORIGINAL_PDB_DIR = env_path("LSP_REFERENCE_DIR", f"{WORKFLOW}/02_rf_diffusion_gai/cyclic_refine")
# 侧链包装PDB的根目录（03_MPNN下的cyclic_youhua_1）
PACKED_ROOT = env_path("LSP_INPUT_DIR", f"{WORKFLOW}/03_MPNN/cyclic_refine")
# 结果保存根目录
DEST_ROOT = env_path("LSP_OUTPUT_DIR", f"{WORKFLOW}/04_select_MPNN/cyclic_refine")
# 日志文件路径
LOG_FILE = f"{DEST_ROOT}/pymol_align.log"

# 关键修复：先创建结果根目录，确保日志文件能正常写入
os.makedirs(DEST_ROOT, exist_ok=True)

# 配置日志（记录成功/失败信息）
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# ============ 2. 初始化PyMOL（静默模式） ============
pymol.finish_launching(['pymol', '-cq'])  # 不显示GUI，安静运行

# ============ 3. 批量处理函数 ============
def process_sample(sample_name, sample_path):
    """处理单个样本：对齐并替换B链"""
    # 样本相关路径
    # 原始PDB文件名应该与样本名的前缀匹配（如cyclic_T100_018_0_g）
    base_sample_name = sample_name.split('_B_')[0]
    original_pdb = os.path.join(ORIGINAL_PDB_DIR, f"{base_sample_name}.pdb")  # 原始PDB（含B链）
    packed_dir = os.path.join(sample_path, "packed")         # 侧链包装PDB目录
    merge_dir = os.path.join(DEST_ROOT, sample_name, "merge_pdb")       # 结果保存目录
    os.makedirs(merge_dir, exist_ok=True)                                # 创建结果目录

    # 检查原始PDB是否存在
    if not os.path.exists(original_pdb):
        msg = f"样本 {sample_name} 缺失原始PDB：{original_pdb}"
        print(f"[ERROR] {msg}")
        logging.error(msg)
        return

    # 获取该样本的所有侧链包装PDB
    packed_pdbs = glob.glob(os.path.join(packed_dir, f"*.pdb"))
    if not packed_pdbs:
        msg = f"样本 {sample_name} 未找到侧链包装PDB：{packed_dir}"
        print(f"[WARNING] {msg}")
        logging.warning(msg)
        return

    # 处理每个侧链包装PDB
    for packed_pdb in packed_pdbs:
        packed_basename = os.path.basename(packed_pdb)
        output_pdb = os.path.join(merge_dir, packed_basename)  # 输出文件名与原始packed保持一致

        # 跳过已处理的文件
        if os.path.exists(output_pdb):
            print(f"[SKIP] 已存在对齐结果：{output_pdb}")
            continue

        try:
            # 清空PyMOL环境
            cmd.reinitialize()

            # 加载原始PDB（作为参考，含目标B链）和侧链包装PDB（待对齐）
            cmd.load(original_pdb, "original")    # 原始结构（保留A链等其他链，B链为参考）
            cmd.load(packed_pdb, "packed")        # 侧链包装结构（待对齐的B链）

            # 核心操作：将packed的B链与original的B链对齐
            # 对齐依据：B链的CA原子（确保主链对齐）
            alignment = cmd.align(
                "packed and chain B and name CA",  # 待对齐的链（侧链包装的B链）
                "original and chain B and name CA" # 参考链（原始的B链）
            )
            print(f"[ALIGN] {sample_name} - 对齐得分：{alignment[0]:.2f}")

            # 删除original的B链（准备替换）
            cmd.remove("original and chain B")

            # 合并结构：original的非B链 + packed的B链
            cmd.create("merged", "original or (packed and chain B)")

            # 保存合并结果
            cmd.save(output_pdb, "merged")

            # 记录成功信息
            msg = f"样本 {sample_name} 处理成功：{packed_basename}"
            print(f"[OK] {msg}")
            logging.info(msg)

        except Exception as e:
            # 记录错误信息
            msg = f"样本 {sample_name} 处理失败 {packed_basename}：{str(e)}"
            print(f"[ERROR] {msg}")
            logging.error(msg)

# ============ 4. 主流程：遍历所有样本 ============
if __name__ == "__main__":
    # 获取所有样本目录（如cyclic_T100_018_0_g_B_chain_emref_1）
    sample_dirs = [
        d for d in glob.glob(os.path.join(PACKED_ROOT, "*"))
        if os.path.isdir(d) and "cyclic_" in os.path.basename(d)
    ]

    total = len(sample_dirs)
    print(f"[INFO] 检测到 {total} 个样本目录，开始处理...")
    logging.info(f"开始处理 {total} 个样本")

    for i, sample_dir in enumerate(sample_dirs, 1):
        sample_name = os.path.basename(sample_dir)
        print(f"\n===== 处理样本 {i}/{total}：{sample_name} =====")
        process_sample(sample_name, sample_dir)

    # 退出PyMOL
    cmd.quit()
    print("\n[INFO] 全部处理完成！结果保存在：", DEST_ROOT)
    logging.info("全部处理完成")
