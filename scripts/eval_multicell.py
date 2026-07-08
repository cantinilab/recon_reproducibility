import pandas as pd
import numpy as np
import tqdm
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import hummuspy
import tqdm
import pickle

from scipy.sparse import csr_matrix, csc_matrix

import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve
import argparse

def compute_pr(
    grn: np.array,
    true: np.array,
    do_auc: bool = True,
    doplot: bool = True,):
    """
    compute_pr computes the precision and recall metrics for the given GRN and true matrix.

    Args:
        grn (np.array): The Gene Regulatory Network matrix, where each element represents the strength of the regulatory relationship between genes.
        true (np.array): The ground truth matrix, where each element indicates the presence (1) or absence (0) of a regulatory relationship.
        do_auc (bool, optional): Whether to compute the Area Under the Precision-Recall Curve (AUPRC). Defaults to True.
        doplot (bool, optional): Whether to plot the precision and recall metrics. Defaults to True.

    Raises:
        ValueError: If the shape of the GRN and the true matrix do not match.

    Returns:
        dict: A dictionary containing precision, recall, and random precision metrics.
    """
    if grn.shape != true.shape:
        raise ValueError("The shape of the GRN and the true matrix do not match.")
    metrics = {}
    if isinstance(grn, (csr_matrix, csc_matrix)):
        grn = grn.toarray()
    if isinstance(true, (csr_matrix, csc_matrix)):
        true = true.toarray()
    true = true.astype(bool)
    tot = (grn.shape[0] * grn.shape[1]) - grn.shape[0]
    precision = (grn[true] != 0).sum() / (grn != 0).sum()
    recall = (grn[true] != 0).sum() / true.sum()
    rand_prec = true.sum() / tot

    if doplot:
        print(
            "precision: ",
            precision,
            "\nrecall: ",
            recall,
            "\nrandom precision:",
            rand_prec,
        )
    metrics.update(
        {
            "precision": precision,
            "recall": recall,
            "rand_precision": rand_prec,
        }
    )
    # Initialize lists to store precision and recall values
    precision_list = [precision]
    recall_list = [recall]
    # Define the thresholds to vary
    thresholds = np.append(
        np.linspace(0, 1, 101)[:-2], np.log10(np.logspace(0.99, 1, 30))
    )
    thresholds = np.quantile(grn, thresholds)
    # Calculate precision and recall for each threshold
    if do_auc:
        for threshold in tqdm.tqdm(thresholds[1:]):
            precision = (grn[true] > threshold).sum() / (grn > threshold).sum()
            recall = (grn[true] > threshold).sum() / true.sum()
            precision_list.append(precision)
            recall_list.append(recall)
        # Calculate AUPRC by integrating the precision-recall curve
        if 1.0 not in recall_list:
            precision_list.insert(0, rand_prec)
            recall_list.insert(0, recall_list[0])
            precision_list.insert(0, rand_prec)
            recall_list.insert(0, 1.0)
        precision_list = np.nan_to_num(np.array(precision_list))
        recall_list = np.nan_to_num(np.array(recall_list))
        auprc = -np.trapz(precision_list, recall_list)
        metrics["auprc"] = auprc

        # Compute Average Precision (AP) manually
        sorted_indices = np.argsort(-grn.flatten())
        sorted_true = true.flatten()[sorted_indices]

        tp_cumsum = np.cumsum(sorted_true)
        fp_cumsum = np.cumsum(~sorted_true)

        precision_at_k = tp_cumsum / (tp_cumsum + fp_cumsum)
        recall_at_k = tp_cumsum / true.sum()

        ap = np.sum(precision_at_k[1:] * np.diff(recall_at_k))
        metrics["ap"] = ap
        if doplot:
            print("Average Precision (AP): ", ap)
        if doplot:
            print("Area Under Precision-Recall Curve (AUPRC): ", auprc)

    # compute EPR
    # get the indices of the topK highest values in "grn"
    if isinstance(grn, csr_matrix):
        grn = grn.toarray()
    if isinstance(grn, csc_matrix):
        grn = grn.toarray()
    indices = np.argpartition(
        grn.flatten(),
        -int(true.sum()))[-int(true.sum()):]
    # Compute the odds ratio
    true_positive = true[np.unravel_index(indices, true.shape)].sum()
    false_positive = true.sum() - true_positive
    # this is normal as we compute on the same number of pred_pos as true_pos
    false_negative = true.sum() - true_positive
    true_negative = tot - true_positive - false_positive - false_negative
    # Avoid division by zero
    # this is a debugger line
    if true_negative == 0 or false_positive == 0:
        odds_ratio = float("inf")
    else:
        odds_ratio = (true_positive * true_negative) / (
            false_positive * false_negative)

    metrics.update({"epr": odds_ratio})
    if doplot:
        print("EPR:", odds_ratio)
        plt.figure(figsize=(10, 8))
        plt.plot(
            recall_list,
            precision_list,
            marker=".",
            linestyle="-",
            color="b",
            label="p-r",
        )
        plt.plot(
            [recall_list[0], recall_list[-1]],
            [rand_prec, rand_prec],
            linestyle="--",
            color="r",
            label="Random Precision",
        )
        plt.legend(loc="lower left")
        plt.title("Precision-Recall Curve")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.xscale("log")
        plt.grid(True)
        plt.show()
    return metrics


def compute_spearman(
    predictions, perturbations, pvals_perturbations,
    pval_threshold=0.05
):

    predictions, perturbations, pvals_perturbations =\
        filter_on_genes_celltypes(
            predictions, perturbations, pvals_perturbations)

    logF2 = perturbations.copy()
    logF2.iloc[:, 2:] = np.abs(logF2.iloc[:, 2:])
    pvals_perturbations.iloc[:, 2:] = 1  # pvals_perturbations.iloc[:,2:]<pval_threshold
    logF2.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:].values\
        * logF2.iloc[:, 2:].values

    # spearman per perturbation
    coefs, pvals = [], []
    perturbation_names = [
        col
        for col in predictions.columns
        if col not in ["celltype", "gene"] and col in perturbations]
    for perturbation in perturbation_names:
        y = logF2.set_index(["gene", "celltype"], inplace=False)[perturbation]
        y_pred = predictions.set_index(
            ["gene", "celltype"], inplace=False)[perturbation]

        common_idx = y.index.intersection(y_pred.index)
        y = y[common_idx]
        y_pred = y_pred[common_idx]

        corr_coef = scipy.stats.spearmanr(y, y_pred)
        coefs.append(corr_coef.statistic)
        pvals.append(corr_coef.pvalue)
    return coefs, pvals


def compute_spearman_celltypes(
    predictions, perturbations, pvals_perturbations,
    pval_threshold=0.05
):

    predictions, perturbations, pvals_perturbations = \
        filter_on_genes_celltypes(
            predictions, perturbations, pvals_perturbations)

    logF2 = perturbations.copy()
    logF2.iloc[:, 2:] = np.abs(logF2.iloc[:, 2:])
    pvals_perturbations.iloc[:, 2:] = 1  # pvals_perturbations.iloc[:,2:]<pval_threshold
    logF2.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:].values\
        * logF2.iloc[:, 2:].values

    # spearman per perturbation
    coefs, pvals = [], []
    perturbation_names = [
        col
        for col in predictions.columns
        if col not in ["celltype", "gene"] and col in perturbations]
    for perturbation in perturbation_names:
        for celltype in celltypes:
            y = logF2[logF2["celltype"] == celltype].set_index(
                ["gene", "celltype"], inplace=False)[perturbation]
            if (y.values == 0).all():
                continue

            y_pred = predictions[
                predictions["celltype"] == celltype].set_index(
                ["gene", "celltype"], inplace=False)[perturbation]

            common_idx = y.index.intersection(y_pred.index)
            y = y[common_idx]
            y_pred = y_pred[common_idx]

            corr_coef = scipy.stats.spearmanr(y, y_pred)
            coefs.append(corr_coef.statistic)
            pvals.append(corr_coef.pvalue)
    return coefs, pvals


def compute_pr_perturbations(
    predictions, perturbations, pvals_perturbations,
    display=True,
    pval_threshold=0.01, logfc_threshold=1,
    doplot=False
):

    predictions, perturbations, pvals_perturbations = \
        filter_on_genes_celltypes(
            predictions, perturbations, pvals_perturbations)

    logF2 = perturbations.copy()
    logF2.iloc[:, 2:] = np.abs(logF2.iloc[:, 2:])>logfc_threshold
    pvals_perturbations.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:]\
        < pval_threshold
    logF2.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:].values\
        * logF2.iloc[:, 2:].values

    # 
    perturbation_names = [
        col
        for col in predictions.columns
        if col not in ["celltype", "gene"] and col in perturbations]
    metrics = {}
    for perturbation in perturbation_names:
        y = logF2.set_index(
            ["gene", "celltype"], inplace=False)[[perturbation]]
        y_pred = predictions.set_index(
            ["gene", "celltype"], inplace=False)[[perturbation]]

        common_idx = y.index.intersection(y_pred.index)
        y = y.loc[common_idx, :]
        y_pred = y_pred.loc[common_idx, :]

        metrics[perturbation] = compute_pr(
            y_pred.values.T, y.values.T, doplot=doplot)
    return metrics


def compute_pr_perturbations_celltypes(
    predictions, perturbations, pvals_perturbations,
    display=True,
    pval_threshold=0.01, logfc_threshold=1,
    doplot=False
):

    predictions, perturbations, pvals_perturbations = \
        filter_on_genes_celltypes(
            predictions, perturbations, pvals_perturbations)

    logF2 = perturbations.copy()
    logF2.iloc[:, 2:] = np.abs(logF2.iloc[:, 2:])>logfc_threshold
    pvals_perturbations.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:]\
        < pval_threshold
    logF2.iloc[:, 2:] = pvals_perturbations.iloc[:, 2:].values\
        * logF2.iloc[:, 2:].values

    #
    perturbation_names = [
        col
        for col in predictions.columns
        if col not in ["celltype", "gene"] and col in perturbations]
    metrics = {}
    for perturbation in perturbation_names:
        for celltype in logF2["celltype"].unique():
            y = logF2.loc[logF2["celltype"] == celltype, :].set_index(
                ["gene", "celltype"], inplace=False)[[perturbation]]
            y_pred = predictions.loc[
                predictions["celltype"] == celltype, :].set_index(
                    ["gene", "celltype"], inplace=False)[[perturbation]]

            common_idx = y.index.intersection(y_pred.index)
            y = y.loc[common_idx, :]
            y_pred = y_pred.loc[common_idx, :]

            metrics[perturbation + "_" + celltype] = compute_pr(
                y_pred.values.T, y.values.T, doplot=doplot)
    return metrics


def filter_on_genes_celltypes(
    predictions,
    perturbations,
    pvals_perturbations
):
    predictions = predictions.copy()
    perturbations = perturbations.copy()
    pvals_perturbations = pvals_perturbations.copy()
    genes = np.intersect1d(
        predictions["gene"].unique(), perturbations["gene"].unique())
    predictions = predictions[predictions["gene"].isin(genes)]
    perturbations = perturbations[perturbations["gene"].isin(genes)]
    pvals_perturbations = pvals_perturbations[pvals_perturbations["gene"].isin(
        genes)]
    perturbations = perturbations.sort_values(
        by=["gene", "celltype"]).reset_index(drop=True)
    pvals_perturbations = pvals_perturbations.sort_values(
        by=["gene", "celltype"]).reset_index(drop=True)
    predictions = predictions.sort_values(
        by=["gene", "celltype"]).reset_index(drop=True)

    return predictions, perturbations, pvals_perturbations


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 'True', 'TRUE'):
        return True
    elif v.lower() in ('no', 'false', 'False', 'FALSE'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
