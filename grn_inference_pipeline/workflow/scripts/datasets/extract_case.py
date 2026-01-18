import numpy as np
import pandas as pd
import snapatac2 as snap
import scanpy as sc
import mudata as md
import scanpy.external as sce
from scipy.sparse import issparse
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--celltypes', required=True)
parser.add_argument('-s','--n_sample', required=True)
parser.add_argument('-d','--seed', required=True)
parser.add_argument('-g','--n_hvg', required=True)
parser.add_argument('-r','--n_hvr', required=True)
parser.add_argument('-t','--root', required=True)
parser.add_argument('-p','--paired', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
celltypes = args['celltypes']
n_sample = int(args['n_sample'])
seed = int(args['seed'])
n_hvg = int(args['n_hvg'])
n_hvr = int(args['n_hvr'])
root = args['root']
paired = args['paired'].lower() == 'true'
path_output = args['path_output']

print(f"Dataset is {'PAIRED' if paired else 'UNPAIRED'}")

# Read
mdata = md.read_h5mu(path_input)

# Filter celltypes - handle both paired and unpaired
if celltypes != 'all':
    celltypes = celltypes.split(';')
    if paired:
        # Paired data - filter mdata shared obs
        mdata = mdata[np.isin(mdata.obs['celltype'], celltypes)].copy()
        mdata.obs['celltype'] = mdata.obs['celltype'].cat.remove_unused_categories()
    else:
        # Unpaired data - filter each modality
        if 'celltype' in mdata.mod['rna'].obs.columns:
            rna_mask = np.isin(mdata.mod['rna'].obs['celltype'], celltypes)
            mdata.mod['rna'] = mdata.mod['rna'][rna_mask].copy()
            mdata.mod['rna'].obs['celltype'] = mdata.mod['rna'].obs['celltype'].cat.remove_unused_categories()
        if 'celltype' in mdata.mod['atac'].obs.columns:
            atac_mask = np.isin(mdata.mod['atac'].obs['celltype'], celltypes)
            mdata.mod['atac'] = mdata.mod['atac'][atac_mask].copy()
            mdata.mod['atac'].obs['celltype'] = mdata.mod['atac'].obs['celltype'].cat.remove_unused_categories()

# if not celltype in mdata.obs (for paired data)
if paired and 'celltype' not in mdata.obs.columns:
    mdata.obs['celltype'] = 'unknown'
    mdata.obs['celltype'] = mdata.obs['celltype'].astype('category')

# Downsample
if n_sample > 0:
    if paired:
        # Paired data - downsample shared obs
        n_sample = np.min([n_sample, mdata.obs.shape[0]])
        barcodes = mdata.obs.sample(n=n_sample, random_state=seed, replace=False).index
        mdata = mdata[barcodes, :].copy()
    else:
        # Unpaired data - downsample RNA (smaller modality)
        n_sample_rna = np.min([n_sample, mdata.mod['rna'].obs.shape[0]])
        barcodes_rna = mdata.mod['rna'].obs.sample(n=n_sample_rna, random_state=seed, replace=False).index
        mdata.mod['rna'] = mdata.mod['rna'][barcodes_rna, :].copy()

# Extract
rna = mdata.mod['rna']
atac = mdata.mod['atac']
atac.var['chr_'] = [b.split('-')[0] for b in atac.var_names]
atac.var_names = atac.var_names.values


#check sparse
if not issparse(rna.X):
    raise ValueError("RNA data matrix is not sparse.")
if not issparse(atac.X):
    raise ValueError("ATAC data matrix is not sparse.")

# Make sure enough features
rna = rna[:, np.sum(rna.X.toarray() != 0., axis=0) > 3].copy()
atac = atac[:, np.sum(atac.X.toarray() != 0., axis=0) > 3].copy()

print(f"After feature filtering - RNA: {rna.shape}, ATAC: {atac.shape}")

if atac.shape[1] == 0:
    raise ValueError(f"ATAC has no features after filtering. All features had <=3 non-zero counts.")
if atac.shape[0] == 0:
    raise ValueError(f"ATAC has no cells after filtering.")

# Normalize
rna.layers['counts'] = rna.X.copy()
atac.layers['counts'] = atac.X.copy()
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)

# HVG
def filter_hvg(adata, n_hvg):
    sc.pp.highly_variable_genes(adata, batch_key='batch')
    hvg = adata.var.sort_values('highly_variable_nbatches', ascending=False).head(n_hvg).index
    del adata.var['highly_variable']
    del adata.var['means']
    del adata.var['dispersions']
    del adata.var['dispersions_norm']
    del adata.var['highly_variable_nbatches']
    del adata.var['highly_variable_intersection']
    del adata.uns
    del adata.obs['batch']
    return hvg.values.astype('U')

# Assign batch for HVG calculation
if paired:
    # Paired data - use shared batch
    rna.obs['batch'] = mdata.obs['batch'] if 'batch' in mdata.obs.columns else 'batch'
    atac.obs['batch'] = mdata.obs['batch'] if 'batch' in mdata.obs.columns else 'batch'
else:
    # Unpaired data - use modality-specific batch
    rna.obs['batch'] = rna.obs['batch_rna'] if 'batch_rna' in rna.obs.columns else 'batch'
    atac.obs['batch'] = atac.obs['batch_atac'] if 'batch_atac' in atac.obs.columns else 'batch'

print(f"RNA unique batches: {rna.obs['batch'].nunique()}, ATAC unique batches: {atac.obs['batch'].nunique()}")

hvg = filter_hvg(rna, n_hvg)
hvr = filter_hvg(atac, n_hvr)
rna = rna[:, np.isin(rna.var_names.values.astype('U'), hvg)].copy()
atac = atac[:, np.isin(atac.var_names.values.astype('U'), hvr)].copy()

# Filter cells
rna = rna[(rna.X.A != 0).sum(1) > 3, :].copy()
atac = atac[(atac.X.A != 0).sum(1) > 3, :].copy()

if paired:
    # PAIRED DATA: Intersect and integrate
    obs_inter = atac.obs_names.intersection(rna.obs_names)
    print(f"Paired data - intersection size: {len(obs_inter)}")
    rna = rna[obs_inter].copy()
    atac = atac[obs_inter].copy()
    
    # Update mdata
    mdata = mdata[obs_inter, :].copy()
    mdata.mod['rna'] = rna
    mdata.mod['atac'] = atac
    
    # Infer latent space
    mdata.obsm['X_spectral'] = snap.tl.multi_spectral([rna, atac], features=None)[1]
    
    # Integrate batches if multiple batches exist
    n_samples = mdata.obs['batch'].unique().size
    if n_samples > 1:
        sce.pp.harmony_integrate(
            mdata,
            key='batch',
            basis='X_spectral',
            adjusted_basis='X_spectral',
            max_iter_harmony=30
        )
    
    # Umap
    sc.pp.neighbors(mdata, use_rep="X_spectral")
    sc.tl.umap(mdata)
    
    # Clean
    del mdata.obsp
else:
    # UNPAIRED DATA: Keep separate, no intersection, no integration
    print(f"Unpaired data - RNA: {rna.shape[0]} cells, ATAC: {atac.shape[0]} cells")
    mdata.mod['rna'] = rna
    mdata.mod['atac'] = atac
    # No spectral integration, no batch integration, no UMAP for unpaired data

# Desparsify
if issparse(rna.X):
    rna.X = rna.X.A
if issparse(atac.X):
    atac.X = atac.X.A

# Update mdata
mdata.mod['rna'] = rna
mdata.mod['atac'] = atac
mdata.update()

if root != 'None':
    sc.pp.neighbors(mdata, use_rep='X_spectral')
    sc.tl.paga(mdata, groups='celltype')
    mdata.uns['iroot'] = np.flatnonzero(mdata.obs['celltype']  == root)[0]
    sc.tl.dpt(mdata)
    mdata.uns = dict()
    del mdata.obsm['X_diffmap']
    del mdata.obsp

# Save
mdata.write(path_output)
