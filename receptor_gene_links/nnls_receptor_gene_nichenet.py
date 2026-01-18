import pandas as pd
import numpy as np
import scipy.optimize
import tqdm
from joblib import Parallel, delayed

    # Define the function to run nnls for a single column
def nnls_column(i):
    return scipy.optimize.nnls(LR, LG[:, i])[0]


njobs = 75

# loop on organisms
for organism in ["mouse", "human"]:
    ligand_receptor_path = "lr_network_" + organism + "_21122021.tsv"
    ligand_target_path = "ligand_target_matrix_nsga2r_final_" + organism + ".tsv"
    # ligand x gene matrix
    ligand_target = pd.read_csv(ligand_target_path, sep=' ').transpose()
    # ligand, receptor - 2 columns dataframe
    ligand_receptor = pd.read_csv(ligand_receptor_path, sep=' ')
    ligand_receptor = ligand_receptor.drop_duplicates(subset=["from", "to"], inplace=False)
    print(ligand_receptor)
    print(ligand_target.index)
    # one-hot encode as a matrix the ligand_receptor dataframe
    ligand_receptor["score"] = 1
    ligand_receptor = ligand_receptor.pivot(index="to", columns="from", values="score")
    ligand_receptor = ligand_receptor.fillna(0)
    print(ligand_receptor.shape)
    common_ligands = list(set.intersection(set(ligand_receptor.columns), set(ligand_target.index)))
    ligand_receptor = ligand_receptor.loc[:, common_ligands]
    ligand_target = ligand_target.loc[common_ligands, :]


    LR = ligand_receptor.transpose().values
    LG = ligand_target.values
    # Parallel execution with joblib
    RG = Parallel(n_jobs=njobs)(delayed(nnls_column)(i) for i in tqdm.tqdm(range(LG.shape[1])))
    # Convert the result to a matrix
    RG = np.column_stack(RG)

    receptor_gene = pd.DataFrame(RG, index=ligand_receptor.index, columns=ligand_target.columns)

    # multiply receptor-ligand * ligand-gene
#    receptor_gene = ligand_receptor.loc[:, common_ligands] @ ligand_target.loc[common_ligands, :]

    receptor_gene = receptor_gene.stack().reset_index()
    receptor_gene.columns = ['receptor', 'gene', 'score']

    # Filter out zero scores
    receptor_gene = receptor_gene[receptor_gene['score'] != 0]

    receptor_gene.to_csv("nnls_receptor_gene_" + organism + ".tsv", sep='\t')
