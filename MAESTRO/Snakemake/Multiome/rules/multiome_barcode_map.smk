# rule for mapping barcodes between scRNA and scATAC


rule multiome_barcode_map:
    input:
        atac_barcode = config["barcode"]["atac_whitelist"],
        rna_barcode = config["barcode"]["rna_whitelist"],
        atac_filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
        rna_filtercount = "Result/RNA/QC/%s_filtered_gene_count.h5" %(config["outprefix"])
    output:
        all_jointfilter_count = "Result/Multiome/%s_multiome_feature_count.h5" %(config["outprefix"]),
        peak_jointfilter_count = "Result/Multiome/%s_multiome_peak_count.h5" %(config["outprefix"]),
        gene_jointfilter_count = "Result/Multiome/%s_multiome_gene_count.h5" %(config["outprefix"])
    params:
        species = config["species"],
        outdir = "Result/Multiome",
        outpre = config["outprefix"]
    benchmark:
        "Result/Benchmark/%s_multiome_BarcodeMap.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO multiome-barcode-map --format h5 --peakcount {input.atac_filtercount} --genecount {input.rna_filtercount} "
        "--atac-whitelist {input.atac_barcode} --rna-whitelist {input.rna_barcode} --species {params.species} "
        "-d {params.outdir} --outprefix {params.outpre}"
