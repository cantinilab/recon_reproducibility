import matplotlib
import pandas as pd
import muon as mu
import scanpy as sc
import numpy as np
import argparse
import decoupler as dc

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-r', '--rna', type=str, required=True)
parser.add_argument('-a', '--atac', type=str, required=True)
parser.add_argument('-m', '--mdata', type=str, required=True)
parser.add_argument('-c', '--n_cores', type=int, required=True)
args = parser.parse_args()

# Load the data
rna = sc.read_10x_h5(args.rna)
rna.var_names_make_unique()
atac = sc.read_h5ad(args.atac)


######################
# Preprocess RNA-seq #
######################

# Basic filtering
sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)

# Annotate the group of mitochondrial genes as 'mt'
rna.var['mt'] = rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

## Filter cells following standard QC criteria.
rna = rna[rna.obs.n_genes_by_counts < 2500, :]
rna = rna[rna.obs.pct_counts_mt < 5, :]

# Normalize the data
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
rna.layers['log_norm'] = rna.X.copy()

sc.pp.highly_variable_genes(rna, min_mean=0.0125, max_mean=3, min_disp=0.5)
# Regress and scale the data
sc.pp.regress_out(rna, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack')
# Restore X to be norm counts
dc.swap_layer(rna, 'log_norm', X_layer_key=None, inplace=True)
# Compute distances in the PCA space, and find cell neighbors
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=40)
# Generate UMAP features
sc.tl.umap(rna)
# Run leiden clustering algorithm
sc.tl.leiden(rna)

# Get cell type markers
markers = dc.get_resource('PanglaoDB')  # Query Omnipath and get PanglaoDB
# Filter by canonical_marker and mouse
markers = markers[markers["mouse"].astype(bool) & markers['canonical_marker'] & (markers['mouse_sensitivity'].astype(float) > 0.5) & (markers['organ'] == "Immune system")]
markers = markers[~markers.duplicated(['cell_type', 'genesymbol'])]
markers['genesymbol'] = markers['genesymbol'].str.capitalize()

print(markers)
print(rna.var_names)

dc.run_ora(
    mat=rna,
    net=markers,
    source='cell_type',
    target='genesymbol',
    min_n=3,
    verbose=True,
    use_raw=False,
#    n_jobs=args.n_cores
)
acts = dc.get_acts(
    rna,
    obsm_key='ora_estimate')
# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e
df = dc.rank_sources_groups(
    acts,
    groupby='leiden',
    reference='rest',
    method='t-test_overestim_var'
)
annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

# Add cell type column based on annotation
rna.obs['celltype'] = [annotation_dict[clust] for clust in rna.obs['leiden']]
rna.obs["batch"] = "Lymph_node2"

#######################
# Preprocess ATAC-seq #
#######################

# Annotate the ATAC-seq data
atac.obs['celltype'] = atac.obs["cluster_annotation_coarse"]
atac = atac[atac.obs['batch'] == "MD-scATAC_72_Spl"]
#######################

# Save the annotated data
mdata = mu.MuData(
    {"rna": rna,
     "atac": atac}
     )

mdata.write(args.mdata)