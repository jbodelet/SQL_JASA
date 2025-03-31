import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import numpy as np

# Load latent space data
df = pd.read_csv("cancerDataResults.csv")
df['pheno'] = df['pheno'].astype('category')
df['pheno_code'] = df['pheno'].cat.codes

# Load functional clusters
with open("clusteredFunctions_10.json", "r") as file:
    functions, clusters_index, grid_z = json.load(file).values()

df_clust = pd.DataFrame(functions).T
df_clust.index = clusters_index

# KDE histogram drawer
def draw_hist(ax, df, col, **kwargs):
    for cat in df['pheno'].cat.categories:
        sns.kdeplot(df[df['pheno'] == cat][col], label=cat, shade=True, ax=ax, **kwargs)

# Only use the first 10 clusters
clusters = df_clust.index.unique()[:10]

# 5 rows of clusters (2 per row) + 1 row of histograms (2 histograms)
ncols = 2
nrows = 6  # 5 cluster rows + 1 histogram row

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8, 14), squeeze=False,
                         gridspec_kw={'height_ratios': [1]*5 + [0.6]})

axes = axes.flatten()

# Plot 10 clusters (first 5 rows = rows 0–4)
for i, cluster in enumerate(clusters):
    ax = axes[i]
    cluster_data = df_clust.loc[cluster]
    if isinstance(cluster_data, pd.Series):
        cluster_data = cluster_data.to_frame().T

    # Plot each time series
    for _, series in cluster_data.iterrows():
        ax.plot(grid_z, series, color="grey", alpha=0.5)

    ax.plot(grid_z, cluster_data.mean(), color='black', linewidth=2)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_yticks([0])
    ax.set_ylim([-3, 3])
    ax.grid(True)
    ax.tick_params(axis='both', which='both',
                   bottom=False, labelbottom=False,
                   left=True, labelleft=True)

    # Add cluster label
    textstr = f'Cluster {cluster}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.15, 0.85, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='center', horizontalalignment='center', bbox=props)

# Draw 2 histograms in the **bottom row** (axes[10] and axes[11])
for i in [10, 11]:
    draw_hist(axes[i], df, "sql.1", bw_adjust=0.5)
    axes[i].invert_yaxis()
    axes[i].set_xticks([])
    axes[i].set_yticks([])
    axes[i].set_xlabel("")
    axes[i].set_ylabel("")
    axes[i].spines['right'].set_visible(False)
    axes[i].spines['bottom'].set_visible(False)
    axes[i].spines['left'].set_visible(False)

# Show x-axis labels only on histograms
axes[10].tick_params(bottom=True, labelbottom=True)
axes[11].tick_params(bottom=True, labelbottom=True)

# Adjust layout and spacing
fig.tight_layout()
fig.subplots_adjust(hspace=0.1)

# Save and show
fig.savefig("plot_functionalClusters_10.pdf", bbox_inches='tight')
