# rule for mapping barcodes between scRNA and scATAC


rule multiome_barcode_map:
    input:
        atac_barcode = config["barcode_atac"],
        rna_barcode = config["barcode_rna"],
        atac_filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
        rna_filtercount = "Result/RNA/QC/%s_filtered_gene_count.h5" %(config["outprefix"])
    output:
        all_jointfilter_count = "Result/Multiome/%s_filtered_feature_count.h5" %(config["outprefix"])
    params:
        species = config["species"],
        outdir = "Result/Multiome/Analysis",
        outpre = config["outprefix"]
    benchmark:
        "Result/Benchmark/%s_multiome_BarcodeMap.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO multiome-barcode-map --format h5 --peakcount {input.atac_filtercount} --genecount {input.rna_filtercount} "
        "--atac-barcode-lib {input.atac_barcode} --rna-barcode-lib {input.rna_barcode} --species {input.species} "
        "-d {params.outdir} --outprefix {params.outpre}"
