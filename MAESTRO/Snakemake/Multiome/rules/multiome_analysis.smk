# rule for multiome analysis


rule multiome_analysis:
    input:
        atacobj = "Result/ATAC/Analysis/%s_scATAC_Object.rds" %(config["outprefix"]),
        rnaobj = "Result/RNA/Analysis/%s_scRNA_Object.rds" %(config["outprefix"]),
    output:
        mergeobject = "Result/Multiome/%s_multiome_Object.rds" %(config["outprefix"]),
        annoplot = "Result/Multiome/%s_annotated_wsnn.png" %(config["outprefix"]),
    params:
        outpre = config["outprefix"],
        outdir = "Result/Multiome",
    benchmark:
        "Result/Benchmark/%s_multiome_Analysis.benchmark" %(config["outprefix"])
    shell:
        "Rscript " + RSCRIPT_PATH + "/Multiome_pipe.R --atacobj {input.atacobj} --rnaobj {input.rnaobj} --prefix {params.outpre} --outdir {params.outdir}"
