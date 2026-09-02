
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import *

def plot_cell_type_dist_plot(full_dict, grn):

    cytokines = list(full_dict.keys())
    
    # All alpha values
    alphas = list(dict.fromkeys(
        alpha
        for cytokine in cytokines
        for alpha in full_dict[cytokine].keys()
    ))
    
    # Union of all cell types across all cytokines and alpha values
    celltypes = list(dict.fromkeys(
        celltype
        for cytokine in cytokines
        for alpha in full_dict[cytokine].keys()
        for celltype in full_dict[cytokine][alpha].keys()
    ))
    
    n_alphas = len(alphas)
    x = np.arange(len(celltypes))
    width = 0.07
    
    fig, ax = plt.subplots(figsize=(18, 8))

    median_values = pd.DataFrame(columns = ["Alpha", "Celltype", "Median values"])
    
    for i, alpha in enumerate(alphas):
    
        offset = (i - (n_alphas - 1) / 2) * width
    
        for j, celltype in enumerate(celltypes):
    
            values = []
    
            # Collect values from all cytokines
            for cytokine in cytokines:
    
                if alpha not in full_dict[cytokine]:
                    continue
    
                if celltype not in full_dict[cytokine][alpha]:
                    continue
    
                data = full_dict[cytokine][alpha][celltype]
    
                if data is None:
                    continue
    
                # Single value
                if np.isscalar(data):
                    if not np.isnan(data):
                        values.append(float(data))
    
                # Multiple values
                else:
                    data = np.asarray(data, dtype=float)
                    data = data[~np.isnan(data)]
                    values.extend(data.tolist())
    
            # Skip empty groups
            if len(values) == 0:
                continue
    
            bp = ax.boxplot(
                values,
                positions=[x[j] + offset],
                widths=width * 0.8,
                patch_artist=True,
                showfliers=False
            )
    
            # Color by alpha
            for box in bp["boxes"]:
                box.set_facecolor(f"C{i}")
    
            # Add legend only once per alpha
            if j == 0:
                bp["boxes"][0].set_label(f"α = {alpha}")

            #print (alpha, "--", grn, "--", np.median(values))
            median_values.loc[len(median_values)] = [alpha, celltype, values]
            
    # Axis labels and title (bold)
    name_mapping = {
        "hummus": "HuMMuS",
        "scenic": "SCENIC",
        "scenic_plus": "SCENIC+",
        "dictys": "Dictys"
    }

    
    ax.set_xlabel("Cell type", fontsize=14, fontweight="bold")
    ax.set_ylabel("Value", fontsize=14, fontweight="bold")
    ax.set_title(
        "Distribution across all celltypes using GRN : " + name_mapping[grn],
        fontsize=16,
        fontweight="bold"
    )
    
    # X-axis labels
    ax.set_xticks(x)
    ax.set_xticklabels(
        celltypes,
        rotation=60,
        ha="right",
        fontsize=12,
        fontweight="bold"
    )
    
    # Y-axis labels
    for label in ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontsize(12)
    
    # Legend
    legend = ax.legend(
        title="Alpha",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=12
    )
    
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_fontsize(12)
    
    for text in legend.get_texts():
        text.set_fontweight("bold")
    
    # Thicker axes
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    # Thicker tick marks
    ax.tick_params(axis="both", width=1.5)
    
    plt.tight_layout()
    plt.savefig(grn + ".png", dpi=300, bbox_inches="tight")
    plt.show()

    #return median_values



def plot_grn_distribution(grn_full_list):
        
    grns = list(grn_full_list.keys())
    
    alphas = list(dict.fromkeys(
        alpha
        for grn in grns
        for cytokine in grn_full_list[grn]
        for alpha in grn_full_list[grn][cytokine]
    ))
    
    n_grns = len(grns)
    n_alphas = len(alphas)
    
    # ---------------------------------------------------------
    # X-axis positions
    # ---------------------------------------------------------
    
    x = np.arange(n_grns)
    
    # Width of individual alpha box
    width = 0.7 / n_alphas
    
    fig, ax = plt.subplots(figsize=(18, 8))
    
    # ---------------------------------------------------------
    # Plot
    # ---------------------------------------------------------

    median_values = pd.DataFrame(columns = ["Alpha", "GRN", "Median values"])
    
    for g, grn in enumerate(grns):
    
        for i, alpha in enumerate(alphas):
    
            # Position alpha boxes around the GRN center
            offset = (
                i - (n_alphas - 1) / 2
            ) * width
    
            values = []
    
            # -------------------------------------------------
            # Collect values across ALL cytokines and
            # ALL cell types for this GRN + alpha
            # -------------------------------------------------
    
            for cytokine in grn_full_list[grn]:
    
                if alpha not in grn_full_list[grn][cytokine]:
                    continue
    
                for celltype in grn_full_list[grn][cytokine][alpha]:
    
                    data = grn_full_list[grn][cytokine][alpha][celltype]
    
                    if data is None:
                        continue
    
                    # Single value
                    if np.isscalar(data):
    
                        if not np.isnan(data):
                            values.append(float(data))
    
                    # Multiple values
                    else:
    
                        data = np.asarray(data, dtype=float)
    
                        data = data[~np.isnan(data)]
    
                        values.extend(data.tolist())
    
            # -------------------------------------------------
            # Skip empty groups
            # -------------------------------------------------
    
            if len(values) == 0:
                continue
    
            # -------------------------------------------------
            # Boxplot
            # -------------------------------------------------
    
            bp = ax.boxplot(
                values,
                positions=[x[g] + offset],
                widths=width * 0.8,
                patch_artist=True,
                showfliers=False
            )
    
            # Color by alpha
            for box in bp["boxes"]:
                box.set_facecolor(f"C{i}")
    
            # Add alpha legend only once
            if g == 0:
                bp["boxes"][0].set_label(
                    f"α = {alpha}"
                )
    
            #print (alpha, "--", grn, "--", np.median(values))
            median_values.loc[len(median_values)] = [alpha, grn, values]
    # ---------------------------------------------------------
    # X-axis = GRNs
    # ---------------------------------------------------------
    
    ax.set_xticks(x)
        
    name_mapping = {
        "hummus": "HuMMuS",
        "scenic": "SCENIC",
        "scenic_plus": "SCENIC+",
        "dictys": "Dictys"
    }
    
    new_names = [name_mapping.get(x, x) for x in grns]
    
    ax.set_xticklabels(
        new_names,
        fontsize=12,
        fontweight="bold"
    )
    
    # ---------------------------------------------------------
    # Labels and title
    # ---------------------------------------------------------
    
    ax.set_xlabel(
        "GRN",
        fontsize=14,
        fontweight="bold"
    )
    
    ax.set_ylabel(
        "Value",
        fontsize=14,
        fontweight="bold"
    )
    
    ax.set_title(
        "Distribution across GRNs",
        fontsize=16,
        fontweight="bold"
    )
    
    # ---------------------------------------------------------
    # Y-axis formatting
    # ---------------------------------------------------------
    
    for label in ax.get_yticklabels():
        label.set_fontweight("bold")
        label.set_fontsize(12)
    
    # ---------------------------------------------------------
    # Legend
    # ---------------------------------------------------------
    
    legend = ax.legend(
        title="Alpha",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=12
    )
    
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_fontsize(12)
    
    for text in legend.get_texts():
        text.set_fontweight("bold")
    
    # ---------------------------------------------------------
    # Axes formatting
    # ---------------------------------------------------------
    
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    
    ax.tick_params(
        axis="both",
        width=1.5
    )
    
    plt.tight_layout()
    plt.savefig("grn_distribution.png", dpi=300, bbox_inches="tight")
    plt.show()

    return median_values
    