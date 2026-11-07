#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：
批量处理 RFdiffusion 生成的环肽骨架：
1. 将每个结构的 A 链替换为参考结构 LspA.pdb 的 A 链
2. 对齐参考后保存为 *_g.pdb
3. 输出目录：02_rf_diffusion_gai/cyclic/
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

# ============ 1. 路径配置 ============
###基础路径
BASE_DIR = f"{WORKFLOW}"
###参考蛋白结构
reference_pdb = os.environ.get("LSP_REFERENCE_PDB") or os.environ.get("LSP_TARGET_PDB")
if not reference_pdb:
    raise SystemExit("Set LSP_TARGET_PDB or LSP_REFERENCE_PDB to the receptor chain A PDB")
REF_PDB = env_path("LSP_REFERENCE_PDB", reference_pdb)
###输入目录
INPUT_DIR = env_path("LSP_INPUT_DIR", f"{WORKFLOW}/01_rf_diffusion/cyclic")
###输出目录
OUTPUT_DIR = env_path("LSP_OUTPUT_DIR", f"{WORKFLOW}/02_rf_diffusion_gai/cyclic")
LOG_FILE = os.path.join(OUTPUT_DIR, "pymol_convert.log")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============ 2. 初始化 PyMOL ============
pymol.finish_launching(['pymol', '-cq'])  # -c 安静模式，-q 不显示 GUI

# ============ 3. 扫描输入文件 ============
pdb_files = sorted(glob.glob(f"{INPUT_DIR}/*.pdb"))
print(f"[INFO] 检测到 {len(pdb_files)} 个输入文件")
if not pdb_files:
    raise SystemExit(f"No PDB inputs: {INPUT_DIR}")

# ============ 4. 主循环 ============
with open(LOG_FILE, "a") as log:
    for pdb_path in pdb_files:
        pdb_name = os.path.basename(pdb_path)
        output_path = os.path.join(OUTPUT_DIR, pdb_name.replace(".pdb", "_g.pdb"))

        # 跳过已存在文件
        if os.path.exists(output_path):
            print(f"[SKIP] 已存在：{output_path}")
            continue

        try:
            # 加载参考与目标
            cmd.reinitialize()  # 清空环境
            cmd.load(REF_PDB, "ref")
            cmd.load(pdb_path, "complex")

            # 对齐 A 链
            cmd.align("complex and chain A", "ref and chain A")

            # 删除 complex 的原 A 链
            cmd.remove("complex and chain A")

            # 将 ref 的 A 链复制到 complex
            cmd.create("complex", "complex or (ref and chain A)")

            # 保存
            cmd.save(output_path, "complex")

            log.write(f"成功 {pdb_name}\n")
            print(f"[OK] {pdb_name}")

        except Exception as e:
            log.write(f"失败 {pdb_name}: {e}\n")
            print(f"[ERROR] {pdb_name}: {e}")

# ============ 5. 清理 ============
cmd.quit()
print("[INFO] 全部处理完成 ✅")
