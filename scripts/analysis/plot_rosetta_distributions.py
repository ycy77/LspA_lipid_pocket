
# Shared path settings; scripts can be invoked from any working directory.
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[2] / "lib"))
from lsp_paths import WORKFLOW, DATA, ANALYSIS, RF_ROOT, env_path
_Path(ANALYSIS).mkdir(parents=True, exist_ok=True)

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# 读取最小化后的数据
results = pd.read_csv(env_path('LSP_INPUT_CSV', f'{ANALYSIS}/pyrosetta_youhua_results_minimized2.csv'))

# 过滤掉失败的计算
valid_results = results[results['dSASA'] != 'N/A']
valid_results = valid_results.astype({'dSASA': float, 'ddG': float, 'BUH': float, 'CYCttE': float})

# 剔除CYCttE大于0和ddG大于0的数据点
filtered_results = valid_results[(valid_results['CYCttE'] <= 0) & (valid_results['ddG'] <= 0)]

print(f"有效数据数量: {len(filtered_results)}")

# 定义指标配置 - 使用参考图片中的配色
metrics = [
    {
        'key': 'ddG',
        'label': 'ddG(kcal/mol)',
        'title': 'Rosetta ddG',
        'color': '#4CAF50',  # 绿色
        'panel': 'a'
    },
    {
        'key': 'dSASA',
        'label': 'dSASA(Å²)',
        'title': 'Interface Area',
        'color': '#FF9800',  # 橙色
        'panel': 'b'
    },
    {
        'key': 'BUH',
        'label': 'BUH',
        'title': 'Buried Unsatisfied Hbonds',
        'color': '#2196F3',  # 蓝色
        'panel': 'c'
    },
    {
        'key': 'CYCttE',
        'label': 'CYCttE',
        'title': 'Distribution of total energy',
        'color': '#E91E63',  # 粉色
        'panel': 'd'
    }
]

# 创建4面板图
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharey=False)
axs = axs.flatten()

# 绘制每个面板
for i, metric in enumerate(metrics):
    ax = axs[i]
    data = filtered_results[metric['key']].values

    if len(data) > 0:
        # 计算核密度估计 (KDE)
        kde = gaussian_kde(data)

        # 创建平滑的x轴数据
        x_range = np.linspace(data.min() * 0.95, data.max() * 1.05, 1000)
        density = kde(x_range)

        # 绘制填充的密度曲线
        ax.fill_between(x_range, density, alpha=0.5, color=metric['color'])
        ax.plot(x_range, density, color=metric['color'], linewidth=2)

        # 添加均值线
        mean_val = data.mean()
        ax.axvline(mean_val, color='black', linestyle='--', linewidth=2)

        # 添加图例
        ax.legend([metric['key']], loc='upper right', frameon=False)

        # 设置标题和标签
        ax.set_title(f"{metric['panel']}", fontsize=16, fontweight='bold', loc='left')
        ax.set_title(metric['title'], fontsize=14, loc='center')
        ax.set_xlabel(metric['label'], fontsize=12)
        ax.set_ylabel('Density', fontsize=12)

        # 添加网格
        ax.grid(True, alpha=0.2)

        # 设置科学计数法
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# 调整布局
plt.tight_layout()

# 保存为PDF格式
output_path = f'{ANALYSIS}/density_minimized2.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
print(f"已保存: {output_path}")

# 保存为PNG格式
output_path_png = f'{ANALYSIS}/density_minimized2.png'
plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
print(f"已保存: {output_path_png}")

plt.close()
print("能量最小化后密度分布图绘制完成！")
