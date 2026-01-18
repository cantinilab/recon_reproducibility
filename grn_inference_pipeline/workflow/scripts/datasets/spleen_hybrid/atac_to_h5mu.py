from scipy.io import mmread
import muon as mu
import scanpy as sc
import anndata as ad
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--matrix', type=str, required=True)
parser.add_argument('-v', '--var', type=str, required=True)
parser.add_argument('-p', '--obs', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()

matrix = mmread(args.matrix)
atac = ad.AnnData(X=mmread(args.matrix).tocsr()).T

atac.var = pd.read_csv(args.var, sep='\t', index_col=0)
atac.obs = pd.read_csv(args.obs, sep='\t', index_col=0)
atac.obs["batch"] = atac.obs["Sample"]
atac.write(args.output)
