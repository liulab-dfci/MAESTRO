rule scatac_report:
    input:
        # bulkqc = "Result/QC/" + config["outprefix"] + "_bam_stat.txt",
        qcfrag = "Result/QC/{sample}/{sample}_scATAC_fragment_size.png",
        qcfrip = "Result/QC/{sample}/{sample}_scATAC_cell_filtering.png",
        count = "Result/QC/{sample}/{sample}_filtered_peak_count.h5",
        clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
        # annotateplot = "Result/Analysis/%s_annotated.png" %(config["outprefix"]),
        genescore = "Result/Analysis/{sample}/{sample}_gene_score.h5",
        tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt",
        # acannotateplot = "Result/Analysis/%s_CistromeTop_annotated.png" %(config["outprefix"]),
        # ms4a1trackplot = "Result/Analysis/%s_MS4A1_genetrack.png" %(config["outprefix"]),
        # cd3dtrackplot = "Result/Analysis/%s_CD3D_genetrack.png" %(config["outprefix"]),
    output:
        summaryreport = "Result/Report/{sample}_scATAC_report.html",
    params:
        outpre = "{sample}",
        inputpath = config["input_path"],
        species = config["species"],
        platform = config["platform"],
        inputformat = config["format"],
        outdir = "Result/Report/{sample}",
    benchmark:
        "Result/Benchmark/{sample}_Report.benchmark"
    shell:
        # "cp {input.readdistr} {input.qcmap} {input.qcfrag} {input.qcfrip} {input.clusterplot} {input.annotateplot} {output.outdir};"
        "python " + SCRIPT_PATH + "/scATAC_HTMLReport.py --directory {params.outdir} --outprefix {params.outpre} "
        "--input-path {params.inputpath} --species {params.species} --platform {params.platform} --input-format {params.inputformat}"