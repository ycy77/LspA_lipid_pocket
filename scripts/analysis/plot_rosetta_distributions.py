
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

# Read the minimized-structure data.
results = pd.read_csv(env_path('LSP_INPUT_CSV', f'{ANALYSIS}/pyrosetta_youhua_results_minimized2.csv'))

# Exclude failed calculations.
valid_results = results[results['dSASA'] != 'N/A']
valid_results = valid_results.astype({'dSASA': float, 'ddG': float, 'BUH': float, 'CYCttE': float})

# Exclude data points with CYCttE > 0 or ddG > 0.
filtered_results = valid_results[(valid_results['CYCttE'] <= 0) & (valid_results['ddG'] <= 0)]

print(f"Valid data points: {len(filtered_results)}")

# Configure metrics using the reference color palette.
metrics = [
    {
        'key': 'ddG',
        'label': 'ddG(kcal/mol)',
        'title': 'Rosetta ddG',
        'color': '#4CAF50',  # Green.
        'panel': 'a'
    },
    {
        'key': 'dSASA',
        'label': 'dSASA(Å²)',
        'title': 'Interface Area',
        'color': '#FF9800',  # Orange.
        'panel': 'b'
    },
    {
        'key': 'BUH',
        'label': 'BUH',
        'title': 'Buried Unsatisfied Hbonds',
        'color': '#2196F3',  # Blue.
        'panel': 'c'
    },
    {
        'key': 'CYCttE',
        'label': 'CYCttE',
        'title': 'Distribution of total energy',
        'color': '#E91E63',  # Pink.
        'panel': 'd'
    }
]

# Create a four-panel figure.
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharey=False)
axs = axs.flatten()

# Draw each panel.
for i, metric in enumerate(metrics):
    ax = axs[i]
    data = filtered_results[metric['key']].values

    if len(data) > 0:
        # Compute the kernel density estimate (KDE).
        kde = gaussian_kde(data)

        # Create a smooth x-axis grid.
        x_range = np.linspace(data.min() * 0.95, data.max() * 1.05, 1000)
        density = kde(x_range)

        # Draw the filled density curve.
        ax.fill_between(x_range, density, alpha=0.5, color=metric['color'])
        ax.plot(x_range, density, color=metric['color'], linewidth=2)

        # Add the mean line.
        mean_val = data.mean()
        ax.axvline(mean_val, color='black', linestyle='--', linewidth=2)

        # Add the legend.
        ax.legend([metric['key']], loc='upper right', frameon=False)

        # Set titles and labels.
        ax.set_title(f"{metric['panel']}", fontsize=16, fontweight='bold', loc='left')
        ax.set_title(metric['title'], fontsize=14, loc='center')
        ax.set_xlabel(metric['label'], fontsize=12)
        ax.set_ylabel('Density', fontsize=12)

        # Add the grid.
        ax.grid(True, alpha=0.2)

        # Use scientific notation.
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Adjust the layout.
plt.tight_layout()

# Save as PDF.
output_path = f'{ANALYSIS}/density_minimized2.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight', format='pdf')
print(f"Saved: {output_path}")

# Save as PNG.
output_path_png = f'{ANALYSIS}/density_minimized2.png'
plt.savefig(output_path_png, dpi=300, bbox_inches='tight')
print(f"Saved: {output_path_png}")

plt.close()
print("Post-minimization density plots completed!")
