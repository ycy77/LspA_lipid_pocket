# LspA peptide inhibitor design

LspA 环肽抑制剂设计代码，包含 RF 骨架生成、LigandMPNN 序列设计、HADDOCK3 对接、PRODIGY 评分、结构优化和分析绘图。

## 设计参数

| 阶段 | 环肽长度 | CA / frame 噪声 | 扩散步数 T | 骨架数 |
|---|---|---|---|---|
| 首轮设计 | 8–15 aa | 0.1 / 0.1 | 50–200，间隔10 | 每个T生成100个，共1,600个 |
| 二轮设计 | 8–12 aa | 0.2 / 0.2 | 20、30、40、50 | 每个T生成100个，共400个 |

每个骨架使用 LigandMPNN 生成10条序列。参数集中在 `config/design_presets.json`，也可以用命令行参数覆盖。

## 目录

```text
LspA_lipid_pocket/
├── README.md
├── run.sh                       # 统一运行入口
├── requirements-analysis.txt    # 分析绘图依赖
├── config/                      # 参数和路径配置
├── lib/                         # 共用路径处理
└── scripts/
    ├── design/                  # RF骨架生成
    ├── preparation/             # 受体替换、B链提取和合并
    ├── sequence/                # LigandMPNN及序列筛选
    ├── docking/                 # HADDOCK3对接
    ├── scoring/                 # PRODIGY评分
    ├── refinement/              # OpenMM最小化、Rosetta界面指标
    └── analysis/                # 组成、logo、RMSD、口袋及评分分布
```

## 环境与输入

使用 Python >=3.9 和 Bash。运行各阶段前激活对应工具环境：

| 阶段 | 依赖 |
|---|---|
| 骨架生成 | 支持 `inference.cyclic` 的 RFdiffusion/RFpeptide、模型权重、target adjacency、scaffold库 |
| 序列设计 | LigandMPNN、膜环境模型和侧链模型权重 |
| 结构准备 | PyMOL；单独提取链的脚本只需Python标准库 |
| 对接与评分 | HADDOCK3/CNS、命令行PRODIGY |
| 结构优化 | OpenMM/PDBFixer、PyRosetta |
| 分析绘图 | `pip install -r requirements-analysis.txt` |

结构、模型权重和结果数据放在代码目录之外，通过配置指定路径。

```bash
cd LspA_lipid_pocket
# 新机器上创建配置；已有local.env时直接编辑。
cp config/example.env config/local.env
```

编辑 `config/local.env`，至少设置：

```bash
export RFDIFFUSION_DIR=/path/to/RFdiffusion
export LIGANDMPNN_DIR=/path/to/LigandMPNN
export LSP_TARGET_PDB=/path/to/LspA_chain_A.pdb
export LSP_ADJ_FILE=/path/to/LspA_adj.pt
export LSP_SCAFFOLD_DIR=/path/to/ppi_scaffolds
```

`LSP_TARGET_PDB` 应为保留原受体残基编号的 A 链 PDB。可另设 `LSP_REFERENCE_PDB` 作为替换受体时的参考。二轮模板通过 `--target` 指定包含 A/B 链的复合物。请按本机环境填写 `local.env` 中的路径。

默认输出位于 `results/workflow`、`results/docking`、`results/analysis`；可用 `LSP_WORKFLOW_DIR`、`LSP_DATA_DIR`、`LSP_ANALYSIS_DIR` 修改。下列示例使用默认输出目录，并从代码根目录执行。`local.env` 和生成结果已加入 `.gitignore`。

## 首轮设计

```bash
# 查看帮助；预览一个生成命令。
bash run.sh --help
bash run.sh design/generate_backbones.py --t-values 50 --count-per-t 1 --dry-run

# 生成1,600个骨架，长度8–15 aa，噪声0.1。
bash run.sh design/generate_backbones.py --mode cyclic

# 替换受体A链。
bash run.sh preparation/replace_receptor.py

# 每骨架设计10条序列并包装侧链。
bash run.sh sequence/design_sequences.py \
  --input-dir results/workflow/02_rf_diffusion_gai/cyclic \
  --output-dir results/workflow/03_MPNN/cyclic

# 按最高overall_confidence选出每骨架的代表序列。
bash run.sh sequence/select_best_packed.py \
  --input-dir results/workflow/03_MPNN/cyclic \
  --output-dir results/workflow/04_select_MPNN/cyclic

# 提取B链并连续编号。
bash run.sh preparation/extract_chain.py \
  --input results/workflow/04_select_MPNN/cyclic --pattern '*/packed/*.pdb' \
  --output-dir results/workflow/05_peptide_chains/cyclic --chain B --renumber

# 对接。
bash run.sh docking/run_haddock3_batch.sh

# 对接结构评分。
bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/xiao_haddock3_cyclic_batch/best_pdb_all \
  --summary results/docking/xiao_haddock3_cyclic_batch/global_score_summary.tsv \
  --output-dir results/docking/xiao_haddock3_cyclic_batch
```

## 二轮设计

选定首轮复合物作为模板，生成8–12 aa环肽，噪声0.2。

```bash
bash run.sh design/generate_backbones.py --mode refine \
  --target /path/to/selected_complex.pdb

LSP_INPUT_DIR="$PWD/results/workflow/01_rf_diffusion/cyclic_refine" \
LSP_OUTPUT_DIR="$PWD/results/workflow/02_rf_diffusion_gai/cyclic_refine" \
  bash run.sh preparation/replace_receptor.py

bash run.sh sequence/design_sequences.py \
  --input-dir results/workflow/02_rf_diffusion_gai/cyclic_refine \
  --output-dir results/workflow/03_MPNN/cyclic_refine

bash run.sh sequence/select_best_packed.py \
  --input-dir results/workflow/03_MPNN/cyclic_refine \
  --output-dir results/workflow/04_select_MPNN/cyclic_refine

bash run.sh preparation/extract_chain.py \
  --input results/workflow/04_select_MPNN/cyclic_refine --pattern '*/packed/*.pdb' \
  --output-dir results/workflow/05_peptide_chains/cyclic_refine --chain B --renumber

LSP_INPUT_DIR="$PWD/results/workflow/05_peptide_chains/cyclic_refine" \
LSP_OUTPUT_DIR="$PWD/results/docking/cyclic_refine" \
  bash run.sh docking/run_haddock3_batch.sh

bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/cyclic_refine/best_pdb_all \
  --summary results/docking/cyclic_refine/global_score_summary.tsv \
  --output-dir results/docking/cyclic_refine
```

## 可选的结构优化与分组评分

将二轮全部packed B链合回参考复合物，再最小化和按骨架筛选：

```bash
LSP_OUTPUT_DIR="$PWD/results/workflow/merged/cyclic_refine" \
  bash run.sh preparation/merge_packed_peptides.py

bash run.sh refinement/minimize_cyclic_openmm.py \
  --input-dir results/workflow/merged/cyclic_refine \
  --output-dir results/docking/refined_openmm

bash run.sh scoring/prodigy_batch.py \
  --input-dir results/docking/refined_openmm \
  --pattern '*/merge_pdb/minimized/*.pdb' --grouped \
  --output-dir results/docking/refined_scores

LSP_SUMMARY_TSV="$PWD/results/docking/refined_scores/final_summary.tsv" \
LSP_PDB_ROOT="$PWD/results/docking/refined_openmm" \
LSP_MODEL_SUBDIR="merge_pdb/minimized" \
  bash run.sh refinement/rosetta_interface_metrics.py

bash run.sh analysis/plot_rosetta_distributions.py
```

OpenMM输出保留输入相对目录并追加 `minimized/`；输出目录应放在输入目录之外。分组评分生成 `final_all.tsv`、`final_summary.tsv` 和 `best_structures/`。

## 分析脚本

所有分析结果默认写入 `results/analysis`。

| 脚本 | 主要输入变量 |
|---|---|
| `analysis/plot_peptide_composition.py` | `LSP_INPUT_DIR`：PDB目录；更换输入时使用新的 `LSP_ANALYSIS_DIR` |
| `analysis/plot_sequence_logos.py` | `LSP_INPUT_DIR`：PDB目录 |
| `analysis/calculate_peptide_rmsd.py` | `LSP_INITIAL_DIR`：样本/packed/唯一PDB；`LSP_DOCKED_DIR`：`*_best.pdb`目录，使用原样本命名对应关系 |
| `analysis/analyze_binding_pockets.py` | `LSP_INPUT_DIR`：对接复合物；`LSP_SCORES_TSV`：评分表 |
| `analysis/compare_prodigy_by_length.py` | `LSP_INITIAL_TSV`、`LSP_OPTIMIZED_TSV`：同一骨架优化前后的配对评分 |
| `refinement/rosetta_interface_metrics.py` | `LSP_SUMMARY_TSV`、`LSP_PDB_ROOT`、`LSP_MODEL_SUBDIR`；`LSP_WORKERS`控制进程数 |
| `analysis/plot_rosetta_distributions.py` | `LSP_INPUT_CSV`：含dSASA/ddG/BUH/CYCttE列的CSV |

```bash
LSP_INPUT_DIR=/path/to/complexes \
LSP_ANALYSIS_DIR="$PWD/results/composition" \
  bash run.sh analysis/plot_peptide_composition.py
```

## 运行控制

- RF可用 `--t-values`、`--count-per-t`、`--noise` 和 `--contigs` 覆盖参数；`--dry-run`只预览命令。线性对照使用 `--mode linear`，不进入环肽专用对接流程。
- HADDOCK默认每任务4核、并发1个；通过 `LSP_CPUS_PER_TASK` 和 `LSP_MAX_CONCURRENT` 调整。
- `LSP_PREPARE_ONLY=1` 只生成对接配置，建议配合独立的 `LSP_OUTPUT_DIR`；正式计算使用新目录。
- RF跳过已存在的非空 `_0.pdb`；MPNN和HADDOCK跳过已有任务目录。重新计算时指定新的输出目录。
- 运行入口可通过 `LSP_PYTHON` 指定解释器；各阶段使用对应的已激活环境。
