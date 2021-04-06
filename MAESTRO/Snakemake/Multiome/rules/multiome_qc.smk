# rule for mapping barcodes between scRNA and scATAC


rule multiome_qc:
    input:
        atac_filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
        rna_filtercount = "Result/RNA/QC/%s_filtered_gene_count.h5" %(config["outprefix"]),
        rna_qc = "Result/RNA/QC/%s_count_gene_stat.txt" %(config["outprefix"]),
        atac_qc = "Result/ATAC/QC/singlecell.txt",
    output:
        all_jointfilter_count = "Result/Multiome/%s_joint_filtered_multiome_feature_count.h5" %(config["outprefix"]),
        peak_jointfilter_count = "Result/Multiome/%s_joint_filtered_peak_count.h5" %(config["outprefix"]),
        gene_jointfilter_count = "Result/Multiome/%s_joint_filtered_gene_count.h5" %(config["outprefix"])
    params:
        species = config["species"],
        outdir = "Result/Multiome",
        outpre = config["outprefix"],
        rna_qc = "../RNA/QC/%s_count_gene_stat.txt" %(config["outprefix"]),
        atac_qc = "../ATAC/QC/singlecell.txt",
    benchmark:
        "Result/Benchmark/%s_multiome_QC.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO multiome-qc --format h5 --peakcount {input.atac_filtercount} --genecount {input.rna_filtercount} "
        "--atac-qc {params.atac_qc} --rna-qc {params.rna_qc} --species {params.species} "
        "-d {params.outdir} --outprefix {params.outpre}"
