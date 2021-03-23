# rule for scRNA analysis

rule scrna_analysis:
    input:
        expression = "Result/Multiome/%s_multiome_gene_count.h5" %(config["outprefix"]),
    output:
        specificgene = "Result/RNA/Analysis/%s_DiffGenes.tsv" %(config["outprefix"]),
        clusterplot = "Result/RNA/Analysis/%s_cluster.png" %(config["outprefix"]),
        annotateplot = "Result/RNA/Analysis/%s_annotated.png" %(config["outprefix"]),
        tflist = "Result/RNA/Analysis/%s.PredictedTFTop10.txt" %(config["outprefix"]),
        rnaobject = "Result/RNA/Analysis/%s_scRNA_Object.rds" %(config["outprefix"]),
    params:
        expression = "../../Multiome/%s_multiome_gene_count.h5" %(config["outprefix"]),
        species = config["species"],
        outpre = config["outprefix"],
        outdir = "Result/RNA/Analysis",
        lisadir = config["lisadir"],
        signature = config["signature"]
    benchmark:
        "Result/Benchmark/%s_Analysis.benchmark" %(config["outprefix"])
    threads:
        config["cores"]
    shell:
        "python " + SCRIPT_PATH + "/lisa_path.py --species {params.species} --input {params.lisadir}; "
        "Rscript " + RSCRIPT_PATH + "/scRNAseq_pipe.R --expression {params.expression} --species {params.species} "
        "--prefix {params.outpre} --signature {params.signature} "
        "--outdir {params.outdir} --thread {threads}"

