import pickle
import zlib
from operator import attrgetter
from pathlib import Path, PurePath
from typing import Sequence, Type
from scipy import sparse
import os
import argparse
os.environ['NUMBA_CACHE_DIR'] = '/tmp/'

import numpy as np
import pandas as pd
import muon as mu

import circe as ci

from distributed import LocalCluster, Client
from arboreto.algo import grnboost2, _prepare_input
import loompy as lp

import tqdm
# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--path_mudata', required=True)
parser.add_argument('-l', '--path_loom_rna', required=True)
parser.add_argument('-r', '--path_grnboost2', required=True)
parser.add_argument('-a', '--path_circe', required=True)
parser.add_argument('-c', '--n_cores', required=True, type=int)
parser.add_argument('-n', '--tf_names', required=True)
parser.add_argument('-o', '--organism')
parser.add_argument('-s', '--seed', type=int, default=666)
args = parser.parse_args()

if args.seed is not None:
    seed = args.seed
else:
    # Use the default seed of arboreto if not specified
    seed = 666

args = vars(parser.parse_args())


from arboreto.core import (
    EARLY_STOP_WINDOW_LENGTH,
    RF_KWARGS,
    SGBM_KWARGS,
    infer_partial_network,
    target_gene_indices,
    to_tf_matrix,
)

import argparse
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count


ATTRIBUTE_NAME_CELL_IDENTIFIER = "CellID"
ATTRIBUTE_NAME_GENE = "Gene"
ATTRIBUTE_NAME_REGULONS_AUC = "RegulonsAUC"
ATTRIBUTE_NAME_REGULONS = "Regulons"
ATTRIBUTE_NAME_METADATA = "MetaData"


def load_exp_matrix_as_loom(
    fname,
    return_sparse=False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from loom file.

    :param fname: The name of the loom file to load.
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    if return_sparse:
        with lp.connect(fname, mode="r", validate=False) as ds:
            ex_mtx = ds.layers[""].sparse().T.tocsc()
            genes = pd.Series(ds.ra[attribute_name_gene])
            cells = ds.ca[attribute_name_cell_id]
        return ex_mtx, genes, cells

    else:
        with lp.connect(fname, mode="r", validate=False) as ds:
            # The orientation of the loom file is always:
            #   - Columns represent cells or aggregates of cells
            # 	- Rows represent genes
            return pd.DataFrame(
                data=ds[:, :],
                index=ds.ra[attribute_name_gene],
                columns=ds.ca[attribute_name_cell_id],
            ).T


def load_exp_matrix(
    fname: str,
    transpose: bool = False,
    return_sparse: bool = False,
    attribute_name_cell_id: str = ATTRIBUTE_NAME_CELL_IDENTIFIER,
    attribute_name_gene: str = ATTRIBUTE_NAME_GENE,
) -> pd.DataFrame:
    """
    Load expression matrix from disk.

    Supported file formats are CSV, TSV and LOOM.

    :param fname: The name of the file that contains the expression matrix.
    :param transpose: Is the expression matrix stored as (rows = genes x columns = cells)?
    :param return_sparse: Returns a sparse matrix when loading from loom
    :return: A 2-dimensional dataframe (rows = cells x columns = genes).
    """
    extension = PurePath(fname).suffixes
    if ".loom" in extension:
        return load_exp_matrix_as_loom(
            fname, return_sparse, attribute_name_cell_id, attribute_name_gene
        )
    elif ".h5ad" in extension:
        from anndata import read_h5ad
        adata = read_h5ad(filename=fname, backed="r")
        if return_sparse:
            # expr, gene, cell:
            return adata.X, adata.var_names.values, adata.obs_names.values
        else:
            return pd.DataFrame(
                adata.X.value.todense(),
                index=adata.obs_names.values,
                columns=adata.var_names.values,
            )
    else:
        df = pd.read_csv(
            fname, sep=suffixes_to_separator(extension), header=0, index_col=0
        )
        return df.T if transpose else df


def run_infer_partial_network(
    target_gene_index,
    gene_names,
    ex_matrix,
    tf_matrix,
    tf_matrix_gene_names,
    method_params,
    seed,
):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=seed,
    )
    return n


def run_grnboost2_fast(rna_anndata, tf_names, n_cores=1, method_params = ["GBM", SGBM_KWARGS], seed=666):

    ex_matrix, gene_names, cell_names = load_exp_matrix(rna_anndata, return_sparse=True)
    print(ex_matrix)
    #gene_names = ex_matrix.columns.values


    ex_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_names)

    import dask.array as da
    import scipy.sparse

    # Convert your CSR matrix to a sparse Dask array
#    dask_matrix = da.from_array(ex_matrix.toarray(), chunks="auto")

    # Save the Dask array in .zarr format to retain lazy-loading with sparse chunks
#   dask_matrix.to_zarr("matrix_dask.zarr")
#
    # Reload the sparse Dask array
 #   ex_matrix = da.from_zarr("matrix_dask.zarr")

    with Pool(n_cores) as p:
        adjs = list(
            tqdm.tqdm(
                p.imap(
                    partial(
                        run_infer_partial_network,
                        gene_names=gene_names,
                        ex_matrix=ex_matrix,
                        tf_matrix=tf_matrix,
                        tf_matrix_gene_names=tf_matrix_gene_names,
                        method_params=method_params,
                        seed=seed,
                    ),
                    target_gene_indices(gene_names, target_genes="all"),
                    chunksize=1,
                ),
                total=len(gene_names),
            )
        )
    adj = pd.concat(adjs).sort_values(by="importance", ascending=False)
    return adj



def run_grnboost2(expression_data, tf_names, n_cores=1):

    local_cluster = LocalCluster(n_workers=n_cores, threads_per_worker=1, death_timeout=30)
    custom_client = Client(local_cluster)
#    scattered_data = custom_client.scatter(expression_data)
    rna_network = grnboost2(
        expression_data=expression_data,
        tf_names=tf_names,
        client_or_address=custom_client)
    custom_client.shutdown()
    custom_client.close()

    return rna_network


path_mudata = args['path_mudata']
path_loom_rna = args['path_loom_rna']
path_grnboost2 = args['path_grnboost2']
path_atacnet = args['path_circe']
n_cores = args['n_cores']
tf_names = args['tf_names']
if args['organism'] =='hg38':
    organism = 'human'
elif args['organism'] == 'mm10':
    organism = 'mouse'

# Load the data
mudata = mu.read_h5mu(path_mudata)
rna = mudata["rna"]

# make loom file
row_attrs = { 
    "Gene": np.array(rna.var.index),
}
col_attrs = { 
    "CellID":  np.array(rna.obs.index),
    "nGene": np.array(np.sum(rna.X.transpose() > 0, axis=0)).flatten(),
    "nUMI": np.array(np.sum(rna.X.transpose(), axis=0)).flatten(),
}

lp.create(path_loom_rna,
          rna.X.transpose(),
          row_attrs,
          col_attrs)

if __name__ == '__main__':
    print('Computing networks...')
    tf_names = pd.read_csv(tf_names, header=None)[0].tolist()
    print('Running grnboost2...')
    rna_network = run_grnboost2_fast(
        path_loom_rna,
#        rna_df, 
        tf_names=tf_names,
        n_cores=n_cores,
        method_params = ["GBM", SGBM_KWARGS])
    rna_network.to_csv(path_grnboost2, index=False)
    print('grnboost2 done!')
    print('Running circe')
    atac = ci.add_region_infos(mudata["atac"], sep=('-', '-'))
    ci.compute_atac_network(atac, organism=organism, njobs=n_cores)
    # Create the atacnet network
#    atac = ci.add_region_infos(mudata["atac"], sep=('-', '-'))
#    ci.compute_atac_network(atac, organism=organism, njobs=n_cores)
    print('Circe done!')
    atac_network = ci.extract_atac_links(atac)
    atac_network = atac_network[atac_network["score"]>0]
    atac_network.to_csv(path_atacnet, index=False)

