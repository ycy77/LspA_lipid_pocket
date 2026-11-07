#!/usr/bin/env bash
set -euo pipefail
LSP_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$LSP_ROOT/config/paths.sh"
export LC_ALL=C.UTF-8  # 解决终端乱码

# --------------------------
# 批量环肽对接配置（用户可修改）
# --------------------------
CONDA_ENV="${LSP_HADDOCK_ENV:-}"                             # conda环境（含haddock3）
TARGET_PDB="${LSP_TARGET_PDB:?Set LSP_TARGET_PDB to the receptor chain A PDB}"  # 靶蛋白A链PDB
B_CHAINS_FOLDER="${LSP_INPUT_DIR:-${LSP_WORKFLOW_DIR}/05_peptide_chains/cyclic}"            # 环肽B链目录（所有.pdb放这里）
HOTSPOT_LIST=(129 130 132 136 137 138 140)  # A链热点残基
OUTPUT_ROOT="${LSP_OUTPUT_DIR:-${LSP_DATA_DIR}/xiao_haddock3_cyclic_batch}"       # 批量结果根目录
CPU_PER_TASK="${LSP_CPUS_PER_TASK:-4}"                              # 每任务核数（≥4，避免环肽计算卡顿）
MAX_CONCURRENT="${LSP_MAX_CONCURRENT:-1}"                             # 并发任务数（总占用核数=MAX_CONCURRENT×CPU_PER_TASK）
[[ "$CPU_PER_TASK" =~ ^[1-9][0-9]*$ && "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo 'CPU/concurrency must be positive integers' >&2; exit 2; }
SLEEP_WHEN_FULL=2                          # 并发满时等待秒数
# 环肽专属核心参数（防断键关键）
CLUST_CUTOFF=0.8                            # 聚类阈值（环肽结构相似性）
TOPOAA_CYCLIC_PEPT="true"                   # 强制识别环肽（必须true）
TOPOAA_CYCLIC_DIST="1.8"                    # 肽键检测阈值（1.8Å精准匹配）
TOPOAA_DISULPHIDE_DIST="3.0"                # 二硫键阈值（英式拼写，避免报错）
FLEXREF_W_VDW="0.6"                         # 降低vdW权重（保护环结构）
# --------------------------

# --------------------------
# 初始化变量（任务统计+全局Score文件）
# --------------------------
success_count=0
fail_count=0
total_task=0
GLOBAL_SCORE_FILE="${OUTPUT_ROOT}/global_score_summary.tsv"  # 全局Score统计


# --------------------------
# 环境检查（强化兼容性）
# --------------------------
# Use the active environment, or an explicitly requested Conda environment.
if [[ -n "$CONDA_ENV" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
if [[ "${LSP_PREPARE_ONLY:-0}" != 1 ]]; then
    command -v haddock3 >/dev/null || { echo "Activate a HADDOCK3 environment first" >&2; exit 1; }
fi
# --------------------------
# 基础文件/目录校验（批量专用）
# --------------------------
# 1. 靶蛋白检查
[ -f "${TARGET_PDB}" ] || {
    echo "ERROR 靶蛋白不存在：${TARGET_PDB}"
    exit 1
}

# 2. B链目录+文件数量检查
[ -d "${B_CHAINS_FOLDER}" ] || {
    echo "ERROR B链目录不存在：${B_CHAINS_FOLDER}"
    exit 1
}
shopt -s nullglob
pdb_list=("${B_CHAINS_FOLDER}"/*.pdb)
total_task=${#pdb_list[@]}
shopt -u nullglob
[ "${total_task}" -eq 0 ] && {
    echo "ERROR ${B_CHAINS_FOLDER}目录下无.pdb文件！"
    exit 1
}

# 3. 资源检查（避免过载）
if command -v nproc >/dev/null 2>&1; then
    avail_cpu="$(nproc)"
    need_cpu=$((MAX_CONCURRENT * CPU_PER_TASK))
    if [ "${need_cpu}" -gt "${avail_cpu}" ]; then
        echo "⚠️ 警告：需${need_cpu}核 > 系统${avail_cpu}核，可能卡顿！"
        exit 1
    fi
fi

# 4. 创建结果目录+全局Score文件表头
mkdir -p "${OUTPUT_ROOT}"
real_output_root="$(cd -- "$OUTPUT_ROOT" && pwd)"
[[ -f "$GLOBAL_SCORE_FILE" ]] || printf "Model\tHADDOCK_Score\tStatus\n" > "$GLOBAL_SCORE_FILE"  # 表头

# 日志初始化
echo "========================================"
echo " 环肽批量对接启动（防断键版）"
echo "========================================"
echo " 任务配置：${total_task}个环肽 | 并发${MAX_CONCURRENT}任务 | 每任务${CPU_PER_TASK}核"
echo " 结果根目录：${real_output_root}"
echo " 环肽优化：共价闭环拓扑+首尾残基固定+vdW权重0.6"
echo "========================================"

# --------------------------
# 核心函数：单环肽对接+后处理（防断键+文件整理）
# --------------------------
run_cyclic_batch() {
    local b_pdb_path="$1"
    local b_basename="$(basename -- "$b_pdb_path" .pdb)"  # 配体前缀（如cyclic_T100_002）
    local task_dir="${real_output_root}/${b_basename}"   # 单任务目录
    local log_path="${task_dir}/haddock3.log"
    local status_file="${task_dir}/.status"
    local task_tag="[$b_basename]"
    local b_chain_max_residue  # 环肽最大残基编号
    local emref1_pdb=""        # emref_1.pdb路径
    local emref1_score="NA"    # emref_1的score

    # --------------------------
    # 1. 环肽PDB深度校验
    # --------------------------
    # 检查B链ATOM行数（避免空文件/格式错误）
    local b_atom_count=$(awk '$1 == "ATOM" && $5 == "B" {print $0}' "${b_pdb_path}" | wc -l)
    if [ "${b_atom_count}" -lt 10 ]; then
        echo "${task_tag} ERROR ATOM行数不足（${b_atom_count}行），跳过！"
        echo -e "${b_basename}\t${emref1_score}\tPDB_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    fi

    # 统计环肽残基（排除插入码，确保连续）
    echo "${task_tag}  统计环肽残基..."
    b_chain_max_residue=$(awk '
        $1 == "ATOM" && $5 == "B" {
            clean_resid = $6; gsub(/[^0-9]/, "", clean_resid);
            if (clean_resid ~ /^[0-9]+$/ && clean_resid > 0) print clean_resid
        }
    ' "${b_pdb_path}" | sort -n | tail -n 1)
    # 残基校验（环肽至少3个残基）
    if [ -z "${b_chain_max_residue}" ] || [ "${b_chain_max_residue}" -lt 3 ]; then
        echo "${task_tag} ERROR 残基统计异常（<3个），跳过！"
        echo -e "${b_basename}\t${emref1_score}\tRESID_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    fi
    echo "${task_tag} OK 环肽信息：残基1~${b_chain_max_residue}（共${b_chain_max_residue}个）"

    # --------------------------
    # 2. 清理旧结果+创建目录
    # --------------------------
    if [ -d "${task_dir}" ]; then
        echo "${task_tag} SKIP existing directory (use a new output directory to rerun)."
        return 0
    fi
    mkdir -p "${task_dir}/data" || {
        echo "${task_tag} ERROR 无法创建任务目录"
        echo -e "${b_basename}\t${emref1_score}\tDIR_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }

    # --------------------------
    # 3. 复制核心文件（保护原文件）
    # --------------------------
    cp -- "${b_pdb_path}" "${task_dir}/ligand_B.pdb" || {
        echo "${task_tag} ERROR 复制环肽失败"
        echo -e "${b_basename}\t${emref1_score}\tCOPY_B_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }
    cp -- "${TARGET_PDB}" "${task_dir}/receptor_A.pdb" || {
        echo "${task_tag} ERROR 复制靶蛋白失败"
        echo -e "${b_basename}\t${emref1_score}\tCOPY_TARGET_ERROR" >> "${GLOBAL_SCORE_FILE}"
        fail_count=$((fail_count + 1))
        return 0
    }

    # --------------------------
    # 4. 生成环肽专用约束（闭环+热点）
    # --------------------------
    local ambig_tbl="${task_dir}/ambig.tbl"
    cat > "${ambig_tbl}" <<EOF
# 环肽对接双约束（防断键+定结合）
# 1. 闭环约束：首残N-尾残C（肽键1.36Å±0.05）
assign (resid 1 and segid B and name N) (resid ${b_chain_max_residue} and segid B and name C) 1.36 0.05 0.05
# 2. 热点约束：A链热点-B链重原子（2.0Å±1.0）
EOF
    # 批量添加热点约束
    for resid in "${HOTSPOT_LIST[@]}"; do
        echo "assign (resid ${resid} and segid A and not hydrogen) (segid B and not hydrogen) 2.0 1.0 3.0" >> "${ambig_tbl}"
    done
    echo "${task_tag} OK 生成约束文件：${ambig_tbl}（1闭环+${#HOTSPOT_LIST[@]}热点）"

    # --------------------------
    # 5. 生成环肽专用config.toml（防断键核心配置）
    # --------------------------
    local flex_start=2
    local flex_end=$((b_chain_max_residue - 1))  # 柔性区域：避开首尾残基
    cat > "${task_dir}/config.toml" <<EOF
# HADDOCK3 环肽对接配置（防断键版）
run_dir = "./"
ncores = ${CPU_PER_TASK}
molecules = ["receptor_A.pdb", "ligand_B.pdb"]

# 阶段1：topoaa（生成共价闭环拓扑）
[topoaa]
autohis = false
cyclicpept = ${TOPOAA_CYCLIC_PEPT}          # 全局环肽开关
cyclicpept_dist = ${TOPOAA_CYCLIC_DIST}      # 肽键检测阈值
disulphide_dist = ${TOPOAA_DISULPHIDE_DIST}  # 二硫键阈值（英式拼写）
# A链（非环肽）末端配置
mol1.charged_nter = false
mol1.charged_cter = false
mol1.5_phosphate = false
# B链（环肽）专属配置（核心！）
mol2.charged_nter = false                    # 无游离N端
mol2.charged_cter = false                    # 无游离C端
mol2.5_phosphate = false
mol2.cyclicpept = true                       # 单独指定B链为环肽

# 阶段2：rigidbody（刚性对接，定结合方向）
[rigidbody]
mode = "local"
sampling = 300                             # 足够构象多样性
ambig_fname = "ambig.tbl"
iniseed = 12345                              # 固定种子，可重复
w_vdw = 0.8                                  # 降低排斥力，防环肽被推
w_elec = 0.2
w_air = 2.0                                  # 提高约束权重，定热点

# 阶段3：clustfcc（刚性后聚类，筛优质构象）
[clustfcc]
clust_cutoff = ${CLUST_CUTOFF}
min_population = 5                           # 过滤孤构象
contact_distance_cutoff = 8.0
plot_matrix = false

# 阶段4：seletopclusts（刚性后筛选，减计算量）
[seletopclusts]
top_clusters = 10
top_models = 1
sortby = "score"                             # 按得分选最优

# 阶段5：flexref（柔性对接，保护闭环）
[flexref]
mode = "local"
max_nmodels = 10
ambig_fname = "ambig.tbl"
sampling_factor = 1
tolerance = 1
w_vdw = ${FLEXREF_W_VDW}                     # 0.6防环肽断裂
w_elec = 0.2
w_air = 2.5                                  # 强化约束，保闭环
fle_sta_1 = ${flex_start}                    # 柔性起始：残基2
fle_end_1 = ${flex_end}                      # 柔性结束：残基${b_chain_max_residue}-1
fle_seg_1 = "B"                              # 仅B链柔性

# 阶段6：clustfcc（柔性后聚类，筛稳定构象）
[clustfcc]
clust_cutoff = ${CLUST_CUTOFF}
min_population = 2
contact_distance_cutoff = 8.0
plot_matrix = false

# 阶段7：seletopclusts（柔性后筛选，定最优）
[seletopclusts]
top_clusters = 5
top_models = 1
sortby = "score"

# 阶段8：emref（能量最小化，优化构象）
[emref]
max_nmodels = 5
ambig_fname = "ambig.tbl"
mode = "local"
tolerance = 2
fle_sta_1 = ${flex_start}                    # 柔性区域与flexref一致
fle_end_1 = ${flex_end}
fle_seg_1 = "B"
w_vdw = ${FLEXREF_W_VDW}                     # 保持权重一致，防断键
w_elec = 0.2
w_air = 2.5

# 阶段9：capri（计算指标+提构象，原caprieval修正）
[caprieval]
ncores = ${CPU_PER_TASK}          # 核心数，与任务核数保持一致（如4）
max_cpus = ${CPU_PER_TASK}        # 最大CPU数，与ncores相同
sortby = "score"
sort_ascending = true             # 升序（得分越低越优，符合CAPRI分析逻辑）
clt_threshold = 1           # 默认聚类阈值，无需修改
EOF

    if [[ "${LSP_PREPARE_ONLY:-0}" == 1 ]]; then
        printf 'PREPARED\n' > "$status_file"
        echo "${task_tag} Prepared configuration; no docking launched."
        return 0
    fi

    # --------------------------
    # 6. 运行对接任务
    # --------------------------
    echo "${task_tag}  启动对接（日志：${log_path}）"
    (
        cd "${task_dir}" || exit 1
        export OMP_NUM_THREADS="${CPU_PER_TASK}"  # 绑定核数，防资源竞争
        if haddock3 config.toml --restart 0 > "${log_path}" 2>&1; then
            printf "SUCCESS\n0\n" > "${status_file}"
        else
            printf "FAIL\nhaddock_failed\n" > "${status_file}"
        fi
    )

    # --------------------------
    # 7. 后处理：解压emref_1.pdb + 提score
    # --------------------------
    if [ -f "${status_file}" ] && [ "$(head -n1 "${status_file}")" = "SUCCESS" ]; then
        # 7.1 解压emref_1.pdb.gz到统一目录
        local emref_gz="${task_dir}/7_emref/emref_1.pdb.gz"
        # 统一目录
        BEST_PDB_DIR="${OUTPUT_ROOT}/best_pdb_all"
        mkdir -p "${BEST_PDB_DIR}"

        if [ -f "${emref_gz}" ]; then
            # 解压并重命名为“配体前缀_emref_1.pdb”
            gunzip -c "${emref_gz}" > "${BEST_PDB_DIR}/${b_basename}_emref_1.pdb"
            emref1_pdb="${BEST_PDB_DIR}/${b_basename}_emref_1.pdb"
            echo "${task_tag} OK 解压最优PDB到统一目录：${emref1_pdb}"
        else
            echo "${task_tag} ⚠️ 未找到emref_1.pdb.gz！"
        fi

        # 7.2 从capri_ss.tsv提取emref_1的score）
        local capri_tsv="${task_dir}/8_caprieval/capri_ss.tsv"
        if [ -f "${capri_tsv}" ]; then
            # 用awk匹配emref_1.pdb行，提取score列（第4列）
            emref1_score=$(awk -F'\t' '$1 ~ /emref_1.pdb/ {print $4; exit}' "${capri_tsv}")
            # 校验score是否为数字
            if [[ "${emref1_score}" =~ ^-?[0-9]+\.[0-9]+$ ]]; then
                echo "${task_tag} OK emref_1 score：${emref1_score}"
            else
                emref1_score="NA"
                echo "${task_tag} ⚠️ 提取score失败（非数字）"
            fi
        else
            echo "${task_tag} ⚠️ 未找到capri_ss.tsv！"
        fi

        # 7.3 标记成功+写入全局Score
        success_count=$((success_count + 1))
        echo -e "${b_basename}\t${emref1_score}\tSUCCESS" >> "${GLOBAL_SCORE_FILE}"
        echo "${task_tag}  任务成功！"

    else
        # 任务失败：记录状态
        fail_count=$((fail_count + 1))
        echo -e "${b_basename}\t${emref1_score}\tFAIL" >> "${GLOBAL_SCORE_FILE}"
        echo "${task_tag} ERROR 任务失败（查看日志：${log_path}）"
    fi
}

# --------------------------
# 主循环：批量调度+并发控制
# --------------------------
echo " 开始提交${total_task}个环肽任务..."
pids=()
for bpath in "${pdb_list[@]}"; do
    # 控制并发：当前任务数≥MAX_CONCURRENT时等待
    while true; do
        current_jobs=$(jobs -rp | wc -l)
        if [ "${current_jobs}" -lt "${MAX_CONCURRENT}" ]; then
            break
        fi
        echo " 当前并发：${current_jobs}/${MAX_CONCURRENT}，等待${SLEEP_WHEN_FULL}秒..."
        sleep "${SLEEP_WHEN_FULL}"
    done

    # 启动任务（后台运行，避免阻塞）
    run_cyclic_batch "${bpath}" &
    pids+=("$!")
    :
done

# --------------------------
# 等待所有任务完成+最终统计
# --------------------------
echo "-----------------------------------------"
echo "⏳ 所有任务已提交，等待后台任务完成..."
child_failed=0
for pid in "${pids[@]}"; do
    wait "$pid" || child_failed=1
done

# Background functions cannot update parent counters. Count persisted statuses.
success_count=$(find "$real_output_root" -mindepth 2 -maxdepth 2 -name .status -exec head -n1 {} \; | awk '$0=="SUCCESS" {n++} END {print n+0}')
fail_count=$(awk -F '\t' 'NR>1 && $3!="SUCCESS" {n++} END {print n+0}' "$GLOBAL_SCORE_FILE")
echo "Completed: success=$success_count failed=$fail_count; summary=$GLOBAL_SCORE_FILE"
[[ "$fail_count" -eq 0 && "$child_failed" -eq 0 ]]
