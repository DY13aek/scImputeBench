import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
import os

# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "./data"

# Load benchmark results from CSV
df = pd.read_csv(os.path.join(os.path.dirname(__file__), "../results/benchmark_results.csv"))

# Rank generation
def generate_rank_data(perf_df):
    rank_df = pd.DataFrame()
    rank_df['tool'] = perf_df['tool']
    metrics = [col for col in perf_df.columns if col != 'tool']
    for metric in metrics:
        ascending = 'RMSE' in metric  # logRMSE: lower is better
        rank_df[metric] = perf_df[metric].rank(ascending=ascending, method='min').astype(int)
    return rank_df.to_dict(orient='list')

rank_data = generate_rank_data(df)

import pprint
pprint.pprint(rank_data)

df_perf = df.copy()
df_rank = pd.DataFrame(rank_data)

# Min-Max Scaling
def min_max_scale(series):
    min_val = series.min()
    max_val = series.max()
    if max_val == min_val:
        return pd.Series([0.5] * len(series))
    return (series - min_val) / (max_val - min_val)

# Metric order
metrics = [
    'cell_PCC', 'Pseudobulk_AUPRC', 'AUPRC (cell-level)', 'logRMSE', # Accuracy
    'ARI', 'NMI', 'Purity', 'SCC_FC',                               # Structural
    'Specificity', 'Coexpression', 'DEG F1', 'DEG SCC', 'DEG AUPRC'  # Biological
]

metric_labels = {
    'cell_PCC': 'cell PCC',
    'Pseudobulk_AUPRC': 'Pseudobulk\nAUPRC',
    'AUPRC (cell-level)': 'cell-level\nAUPRC',
    'logRMSE': 'logRMSE',
    'ARI': 'ARI',
    'NMI': 'NMI',
    'Purity': 'Purity',
    'SCC_FC': 'SCC-FC',
    'Specificity': 'Specificity',
    'Coexpression': 'Coexpression',
    'DEG F1': 'F1',
    'DEG SCC': 'SCC',
    'DEG AUPRC': 'AUPRC'
}

category_map = {
    'cell_PCC': 'Accuracy',
    'Pseudobulk_AUPRC': 'Accuracy',
    'AUPRC (cell-level)': 'Accuracy',
    'logRMSE': 'Accuracy',
    'ARI': 'Structural',
    'NMI': 'Structural',
    'Purity': 'Structural',
    'SCC_FC': 'Structural',
    'Specificity': 'Biological',
    'Coexpression': 'Biological',
    'DEG F1': 'Biological',
    'DEG SCC': 'Biological',
    'DEG AUPRC': 'Biological'
}

colors = {
    'Accuracy': '#5e3c99',
    'Structural': '#e66101',
    'Biological': '#2b83ba'
}

# Min-Max Scaling per metric
df_normalized = df_perf.copy()
for metric in metrics:
    df_normalized[metric] = min_max_scale(df_perf[metric])

# Figure
fig, axes = plt.subplots(len(df_perf), len(metrics), figsize=(16, 8),
                         gridspec_kw={'hspace': 0.025, 'wspace': 0.15})

fig.patch.set_alpha(0)

for j, metric in enumerate(metrics):
    for i, tool in enumerate(df_perf['tool']):
        ax = axes[i, j]

        normalized_value = df_normalized.loc[df_normalized['tool'] == tool, metric].values[0]
        rank = df_rank.loc[df_rank['tool'] == tool, metric].values[0]
        rank_normalized = (rank - 1) / 10

        base_color = colors[category_map[metric]]
        rgb = to_rgb(base_color)
        final_color = tuple(1 - rank_normalized * (1 - rgb[k]) for k in range(3))

        ax.barh(0, normalized_value, height=0.7, color=final_color,
                edgecolor='gray', linewidth=0.5)

        if rank <= 3:
            ax.text(-0.08, 0, str(int(rank)), ha='center', va='center',
                   fontsize=13, fontweight='normal', color='black')

        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 0.5)
        ax.set_xticks([])
        ax.set_yticks([])

        for spine in ax.spines.values():
            spine.set_visible(False)

        if j == 0:
            ax.text(-0.35, 0, tool, ha='right', va='center', fontsize=14)

        if i == len(df_perf) - 1:
            ax.text(0.5, -1.2, metric_labels[metric], ha='right', va='top',
                   fontsize=14, rotation=45, transform=ax.transData)

plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.08)

output_path = os.path.join(BASE_DIR, 'benchmark_bargraph_normalized.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.close()

print(f"Bar graph saved to: {output_path}")
print("\nNormalized value range check:")
for metric in metrics:
    print(f"{metric}: {df_normalized[metric].min():.3f} - {df_normalized[metric].max():.3f}")
