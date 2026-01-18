rule download_spleen_hybrid:
    params:
        atac = config['datasets']['spleen_hybrid']['url']['archr_project'],
        rna = config['datasets']['spleen_hybrid']['url']['10X_platform']
    output:
        archr_proj = temp(directory("datasets/spleen_hybrid/atac/archr_mouse_immune_cell_atlas")),
        atac_mtx = temp("datasets/spleen_hybrid/atac/atac.mtx"),
        atac_peaks = temp("datasets/spleen_hybrid/atac/atac_peaks.csv"),
        atac_obs = temp("datasets/spleen_hybrid/atac/atac_obs.csv"), 
        atac = "datasets/spleen_hybrid/atac.h5ad",
        rna = "datasets/spleen_hybrid/rna.h5"
    singularity:
        "workflow/envs/gretabench.sif"
    shell:
        # 1. ATAC-seq: GEO242466
        # 2. RNA-seq: 10x datasets - Mixture of Cells from Mouse Lymph Nodes (Lymph node 1)
        """
        mkdir -p datasets/spleen_hybrid/atac
        mkdir -p datasets/spleen_hybrid/rna
        wget "{params.atac}" -O "{output.atac}.tar.gz"
        tar -xzvf "{output.atac}.tar.gz" -C datasets/spleen_hybrid/atac
        rm "{output.atac}.tar.gz"

        Rscript workflow/scripts/datasets/spleen_hybrid/atac_archr.R \
            {output.archr_proj} \
            {output.atac_mtx} \
            {output.atac_peaks} \
            {output.atac_obs}

        python workflow/scripts/datasets/spleen_hybrid/atac_to_h5mu.py \
            --matrix {output.atac_mtx} \
            --var {output.atac_peaks} \
            --obs {output.atac_obs} \
            --output {output.atac}

        wget {params.rna} -O {output.rna}
        """

rule annotate_spleen_hybrid:
    input:
        atac = "datasets/spleen_hybrid/atac.h5ad",
        rna = "datasets/spleen_hybrid/rna.h5"
    params:
        n_cores = 10
    output:
        out = "datasets/spleen_hybrid/annotate.h5mu"
    singularity:
        "workflow/envs/gretabench.sif"
    shell: # Format as single h5mu file
        # 1. R: Open Seurat atac, save as h5ad
        # 2. python: Open h5ad and others, save as h5mu
        """
        python workflow/scripts/datasets/spleen_hybrid/annotate.py \
            --atac {input.atac} \
            --rna {input.rna} \
            --mdata {output.out} \
            -c {params.n_cores}
        """
