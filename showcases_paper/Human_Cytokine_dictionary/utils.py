
import pandas as pd
import scanpy as sc
import pandas as pd
from recon.explore import Celltype
import numpy as np
import scanpy as sc  # single cell data
import pandas as pd  # data manipulation
import liana as li  # cell communication
import recon  # multilayer and perturbation prediction
import recon.data
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

def bar_plot(aurocs):
    plt.figure(figsize=(7,4))

    bars = plt.bar(aurocs.keys(), aurocs.values())
    
    plt.ylabel("AUROC", fontsize=12, fontweight="bold")
    plt.xlabel("Cell Type", fontsize=12, fontweight="bold")
    plt.ylim(0, 1)
    
    # Bold tick labels
    plt.xticks(rotation=45, ha="right", fontsize=11, fontweight="bold")
    plt.yticks(fontsize=11, fontweight="bold")
    
    # Add AUROC values on top of bars
    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width()/2,
            height + 0.01,
            f"{height:.2f}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold"
        )
    
    plt.tight_layout()
    plt.show()


def alpha_experiment(direct_effect, indirect_effect, cell__types, dct_ct, gt):
        
    dct_alphas = {}
    
    for alpha in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        total_effect = recon.explore.combine_effects(direct_effect, indirect_effect, alpha=alpha)
    
        #aurocs = dict.fromkeys(cell__types, 0)
        aurocs = {}
        
        for ct in cell__types: #dct_ct.keys():
            truth = gt[gt["celltype"] == dct_ct[ct]].copy()
        
            #print (truth.shape)
            pred = total_effect[ct].rename("score").reset_index()
            pred.columns = ["gene", "score"]
        
            #print (pred)
            # Binary labels
            truth["truth"] = (
                (truth["adj_p_value"] < 0.05) &
                (truth["log_fc"] > 0)
            ).astype(int)
        
            truth = truth[["gene", "truth"]]
           
            # If duplicate genes exist, keep the maximum truth value
            truth = truth.groupby("gene", as_index=False)["truth"].max()
            #print (truth) 
            # -------------------------
            # Merge on common genes
            # -------------------------
            
            merged = pred.merge(truth, on="gene", how="inner")
           
            #print(f"Predicted genes : {len(pred)}")
            #print(f"Truth genes     : {len(truth)}")
            #print(f"Common genes    : {len(merged)}")

            if merged.empty:
                #print(f"No common genes")
                auc = np.nan
            else:
                auc = roc_auc_score(
                    merged["truth"],
                    merged["score"]
                )
            
            #print(f"AUROC: {auc:.4f}")
            aurocs[ct] = auc
    
        #bar_plot(aurocs)
        dct_alphas[alpha] = aurocs
        
    return dct_alphas
    
def joint_plot(dct_alphas, title):
    
    df = pd.DataFrame(dct_alphas)
    
    ax = df.plot(
        kind="bar",
        figsize=(18,6),
        width=0.85
    )
    
    ax.set_xlabel("Cell Type", fontsize=14, fontweight="bold")
    ax.set_ylabel("AUROC", fontsize=14, fontweight="bold")
    ax.set_ylim(0, 1)
    
    plt.xticks(rotation=30, ha="right", fontsize=12, fontweight="bold")
    plt.yticks(fontsize=12, fontweight="bold")
    
    plt.legend(
        title="Alpha",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=10,
        title_fontsize=11
    )
    
    # Add values on top of each bar
    for container in ax.containers:
        ax.bar_label(
            container,
            fmt="%.2f",
            fontsize=7,
            fontweight="bold",
            rotation=90,
            padding=2
        )

    plt.title(title + ": AUROC Across Cell Types for Different α Values",
          fontsize=16,
          fontweight="bold")
    
    plt.tight_layout()
    plt.savefig("./Figures/" + title + ".png", dpi=300, bbox_inches="tight")
    plt.show()

    