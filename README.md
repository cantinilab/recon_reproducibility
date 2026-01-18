This repository contains the code to replicate the experiment of ReCoN's manuscript, accessible [here]()

We still need to upload all necessary files for the reproduction on Zenodo, it should be out in the coming days!
<br>Zenodo link: []()

## Organisation

The repository is organised in three subfolders:
### **1. [Notebooks of the showcases](https://github.com/cantinilab/recon_reproducibility/tree/main/showcases_paper)**

The Jupyter Notebook reproduces the main showcases of the paper, the figures, and supplementary figures.
They are organised per dataset, so the [Heart Failure notebook]() contains both the "multicellular coordination" and the "multicellular program drivers" showcases.

### **2. [Receptor gene links](https://github.com/cantinilab/recon_reproducibility/tree/main/receptor_gene_links)**
The code to obtain the receptor-gene links for both human and mouse. It relies on the ligand-receptor and the ligand-target genes networks of [Nichenet]()https://www.nature.com/articles/s41592-019-0667-5, which can be [downloaded on Zenodo]().

### **3. [GRN inference pipeline](https://github.com/cantinilab/recon_reproducibility/tree/main/grn_inference_pipeline)**
A snakemake workflow to reproduce GRNs using HuMMuS, from data downloading to network inference.
It is one of the subcomponents of ReCoN's networks and analyses, but HuMMuS can be tricky to use, so we provide it as a Snakemake pipeline. This pipeline is heavily based on [GRETA](hhtps://github.com/saezlab/GRETA).

⚠️ There is some stochasticity in HuMMuS inference, so the obtained network might slightly differ/ If you need the exact same GRNs than the one used in the publication, please download them [here, on Zenodo]().

⚠️⚠️ It is also a heavy computation step, which might run for more than 10 hours on a classic HPC.

Please refer to the README inside the folder for more details.

## TO DO/missing

- 'multicell' and 'eval_multicell' scripts to be uploaded - in the Notebooks of the showcases folder
- release all Zenodo files - once the preprint is out
- precisions on recon version and environments for Notebooks of the showcases folder
- Fill up links in the README




