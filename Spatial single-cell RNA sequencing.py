import sys
import pandas as pd
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import gc

import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl

# Set paths to data
sp_data_folder = 'PATH/'
sp_data_folder

adata = sc.read_visium(
    sp_data_folder, count_file='/path/', load_images=True)

adata.var_names_make_unique()

adata
print(adata.obs.head())  # We have 2987 such observations (cells)
adata.var.head()  # We have 31053 such variables (genes)


sc.pl.highest_expr_genes(adata, n_top=20, show=False)
plt.title("Top 20 Highly Expressed Genes")
plt.xlabel("Gene Rank")
plt.ylabel("Fraction of Total Counts")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()


#####
if not isinstance(adata.X, np.ndarray):
    X_dense = adata.X.toarray()
else:
    X_dense = adata.X

# Compute mean expression for each gene
gene_means = np.mean(X_dense, axis=0)

# Get top 20 genes by average expression
top_n = 20
top_indices = np.argsort(gene_means)[::-1][:top_n]
top_gene_names = adata.var_names[top_indices]
top_gene_values = gene_means[top_indices]

# Plot
plt.figure(figsize=(12, 6))
plt.bar(top_gene_names, top_gene_values, color='skyblue')
plt.xticks(rotation=90)
plt.title("Top 20 Highly Expressed Genes")
plt.xlabel("Gene")
plt.ylabel("Mean Expression")
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()


####
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.sparse import issparse

def plot_top_expressed_genes(adata, n_genes=20, log_transform=True):
    """
    Plot top expressed genes with optional log2 transformation

    Args:
        adata: AnnData object containing gene expression data
        n_genes: Number of top genes to display
        log_transform: Whether to apply log2 transformation to expression values
    """
    # Convert to dense array if sparse
    X_dense = adata.X.toarray() if issparse(adata.X) else adata.X

    # Apply log2 transformation if requested (adding pseudocount to avoid log(0))
    if log_transform:
        X_dense = np.log2(X_dense + 1)

    # Compute mean expression for each gene
    gene_means = np.mean(X_dense, axis=0)

    # Get top genes by average expression
    top_indices = np.argsort(gene_means)[::-1][:n_genes]
    top_gene_names = adata.var_names[top_indices]
    top_gene_values = gene_means[top_indices]

    # Create plot
    plt.figure(figsize=(12, 6))
    bars = plt.bar(top_gene_names, top_gene_values, color='skyblue')

    # Add value labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                 f'{height:.2f}',
                 ha='center', va='bottom', fontsize=9)

    # Formatting
    plt.xticks(rotation=45, ha='right')
    plt.title(
        f"Top {n_genes} {'Log2 ' if log_transform else ''}Expressed Genes")
    plt.xlabel("Gene")
    plt.ylabel("Mean Expression" + (" (log2)" if log_transform else ""))
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()


# Example usage:
plot_top_expressed_genes(adata, n_genes=20, log_transform=True)

# For log2 transformed values (recommended for visualization)
plot_top_expressed_genes(adata, n_genes=20, log_transform=True)

# For raw counts (without transformation)
plot_top_expressed_genes(adata, n_genes=20, log_transform=False)


# Convert the sparse matrix to a dense matrix
dense_matrix = adata.X.toarray()

# Print the dense matrix
print(dense_matrix)


cell_index = 0  # Index of the cell you want to examine

# Access the expression data for the specific cell
gene_expression = adata.X[cell_index]

# Get the indices of the genes with non-zero expression in the cell
expressed_gene_indices = gene_expression.nonzero()[1]

# Access the gene names from the 'var_names' attribute
expressed_genes = adata.var_names[expressed_gene_indices]

# Print the list of expressed genes in the cell
print(expressed_genes)


# Change alpha to 0 to see the tissue sample or plot by setting color=None
sc.pl.spatial(adata, img_key="hires", alpha=0.5)

####
sc.pl.spatial(
    adata,
    img_key="hires",
    color=["Rorb", "Vip"],   # genes of interest
    cmap="magma",            # better color contrast
    alpha_img=0.6,           # slightly dim background tissue image
    size=1.2,                # size of spatial spots
    vmin=0,
    vmax='p99.2',            # cap extremes for better contrast
    show=True
)
####
Optional Style Tweaks


def plot_spatial_expression(adata, genes, title=None):
    sc.pl.spatial(
        adata,
        img_key="hires",
        color=genes,
        cmap="plasma",
        alpha_img=0.6,
        size=1.2,
        vmin=0,
        vmax='p99.2',
        show=False
    )
    if title:
        plt.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.show()


plot_spatial_expression(adata, ["Vip", "mt-Co3"],
                        title="Spatial Expression of Vip and mt-Co3")
####
sc.pl.spatial(adata, img_key="hires", color=None)

gene_names = adata.var.index

print("Is Rorb gene present in the vars?", "Rorb" in gene_names)
print("Is Vip gene present in the vars?", "Vip" in gene_names)

####

with mpl.rc_context({'figure.figsize': [6, 7],
                     'axes.facecolor': 'black'}):
    sc.pl.spatial(adata, color=["Rorb", "Vip", "mt-Co3"], img_key=None, size=1,
                  vmin=0, cmap='magma', vmax='p99.0',
                  gene_symbols='SYMBOL'
                  )


####

with plt.rc_context({'figure.figsize': (6, 7), 'axes.facecolor': 'black'}):
    # Spatial plot of selected genes with enhanced color and intensity scaling
    sc.pl.spatial(
        adata,
        color=["Rorb", "Vip", "mt-Co3"],  # Genes to visualize
        img_key=None,                     # Set to None if no tissue image is used
        size=1,                           # Dot size
        vmin=0,                           # Minimum value for color scale
        vmax='p99.0',                     # Cap at 99th percentile
        cmap='magma',                     # Colormap for expression
        gene_symbols='SYMBOL'             # Column in adata.var with gene names
    )

####

sc.pp.calculate_qc_metrics(adata, inplace=True)

# adata.var[adata.var.index.str.startswith('mt-')]
adata.var["mt"] = adata.var_names.str.startswith("mt-")

adata.var

ribo_genes = pd.read_csv(
    '/home', skiprows=2, header=None)
ribo_genes

adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)

adata.var

adata

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", 'ribo'], inplace=True)

adata

number_of_spots = adata.obs_names.shape[0]
mean_reads_per_spot = adata.obs['total_counts'].mean()
median_genes_per_spot = adata.obs['n_genes_by_counts'].median()

print("Number of spots under tissue:", number_of_spots)
print("Mean reads per spot:", mean_reads_per_spot)
print("Median genes per spot:", median_genes_per_spot)

adata.var


adata.obs

fig, axs = plt.subplots(1, 2, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])

adata.obs.sort_values('total_counts')

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)


sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20]
adata = adata[adata.obs["pct_counts_ribo"] < 2]

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
             jitter=0.4, multi_panel=True)


# Filter genes based on minimum number of cells.
sc.pp.filter_genes(adata, min_cells=10)
adata.var

# Sum the gene expression counts across cells
gene_counts_sum = np.sum(adata.X, axis=0)
# Find the gene with the highest sum
gene_with_most_counts = adata.var_names[np.argmax(gene_counts_sum)]
print("Gene with the most counts:", gene_with_most_counts)

# Specify the gene name
gene_name = 'mt-Co3'

# Get the index of the gene in adata.var_names
gene_index = list(adata.var_names).index(gene_name)

# Extract the gene expression counts for the gene
gene_counts = adata.X[:, gene_index]

# Convert sparse matrix to array
gene_counts_array = gene_counts.toarray().flatten()

# Plot the histogram
plt.hist(gene_counts_array, bins=30)
plt.title(f"Histogram of Counts per Cell for Gene {gene_name}")
plt.xlabel("Counts per Cell")
plt.ylabel("Frequency")
plt.show()

sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)


# Get the normalized gene expression counts for the gene
gene_counts_normalized = adata.X[:, gene_index]

# Convert sparse matrix to array
gene_counts_normalized_array = gene_counts_normalized.toarray().flatten()

# Plot the histogram
plt.hist(gene_counts_normalized_array, bins=30)
plt.title(f"Histogram of Normalized Counts per Cell for Gene {gene_name}")
plt.xlabel("Normalized Counts per Cell")
plt.ylabel("Frequency")
plt.show()

sc.pp.log1p(adata)


# Get the log-transformed gene expression counts for the gene
gene_counts_log_transformed = adata.X[:, gene_index]

# Convert sparse matrix to array
gene_counts_log_transformed_array = gene_counts_log_transformed.toarray().flatten()

# Plot the histogram
plt.hist(gene_counts_log_transformed_array, bins=30)
plt.title(f"Histogram of Log-Transformed Counts per Cell for Gene {gene_name}")
plt.xlabel("Log-Transformed Counts per Cell")
plt.ylabel("Frequency")
plt.show()

sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

adata.var

sc.pl.highly_variable_genes(adata)


adata = adata[:, adata.var.highly_variable]

sc.pp.pca(adata)  # By default calculates 30 PCAs

sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_pcs=20)


sc.tl.umap(adata, n_components=2)

sc.pl.umap(adata)

sorted_by_num_cells = adata.var['n_cells_by_counts'].sort_values(
    ascending=False)
sorted_by_num_cells


sc.pl.umap(adata, color=["Gm42418", "Fth1", "3110035E14Rik"])

# Experiment with values 0.3, 0.6, and 1.0. View the result in the UMAP plot below.
sc.tl.leiden(adata, resolution=0.6, key_added="clusters")
adata.obs

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts",
           "n_genes_by_counts", "clusters"], wspace=0.4)

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=[
              "total_counts", "n_genes_by_counts"])

####
plt.rcParams.update({
    "figure.figsize": (8, 8),  # Square figure
    "figure.facecolor": "white",  # White background
    "figure.dpi": 100,  # Adjust resolution if needed
})

# Create spatial plot with enhanced parameters
sc.pl.spatial(
    adata,
    img_key="hires",
    color=["total_counts", "n_genes_by_counts"],  # Features to visualize
    ncols=2,  # Display plots in 2 columns
    wspace=0.5,  # Add space between subplots
    frameon=False,  # Cleaner look without frames
    colorbar_loc="right",  # Consistent colorbar placement
    save=None,  # Set to filename if you want to save
    show=True  # Set to False if you're embedding in a larger figure
)
####

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1)


xmin = adata.obsm['spatial'][:, 0].min()
xmax = adata.obsm['spatial'][:, 0].max()
ymin = adata.obsm['spatial'][:, 1].min()
ymax = adata.obsm['spatial'][:, 1].max()

print("x-coordinate range:", xmin, "to", xmax)
print("y-coordinate range:", ymin, "to", ymax)


sc.pl.spatial(adata, img_key="hires", color="clusters", groups=["0", "2", "5"], crop_coord=[
              # coord = [xmin, xmax, ymin, ymax]
              3000, 10000, 5000, 10000], alpha=0.4, size=0.7)
# To Visualize the underlying image set alpha=0

# Group by clusters and perform differential gene expression analysis using t-tets (also try 'wilcoxon')
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups(adata, n_genes=10)

sc.pl.rank_genes_groups_heatmap(
    adata, groups="6", n_genes=10, groupby="clusters")

plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["clusters", "Ctxn1"])

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
top_features = {}
n_top_genes = 10  # desired number of top genes per cluster
for group in groups:
    top_features[group] = result["names"][group][:n_top_genes]

# Print the top features for each cluster
for group, features in top_features.items():
    print(f"Cluster {group} top features:")
    for feature in features:
        print(feature)
    print()


# Access the marker genes results from rank_genes_groups
marker_genes = adata.uns['rank_genes_groups']

# Iterate over each group and print the marker genes
for group in marker_genes['names'].dtype.names:
    print(f"Group: {group}")
    print(marker_genes['names'][group][:10])  # Print the top 10 marker genes
    print("\n")
