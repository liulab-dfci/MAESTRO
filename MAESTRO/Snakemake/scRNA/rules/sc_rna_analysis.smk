
rule scrna_analysis:
    input:
        expression = "Result/QC/{sample}/{sample}_filtered_gene_count.h5"
    output:
        specificgene = "Result/Analysis/{sample}/{sample}_DiffGenes.tsv",
        clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
        annotateplot = "Result/Analysis/{sample}/{sample}_annotated.png",
        tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt"
    params:
        expression = "../QC/{sample}/{sample}_filtered_gene_count.h5",
        species = config["species"],
        outpre = "{sample}/{sample}",
        outdir = "Result/Analysis",
        lisadir = config["lisadir"],
        signature = config["signature"]
    benchmark:
        "Result/Benchmark/{sample}/{sample}_Analysis.benchmark"
    threads:
        config["cores"]
    shell:
        """
        python {SCRIPT_PATH}/lisa_path.py --species {params.species} --input {params.lisadir};
        Rscript {RSCRIPT_PATH}/scRNAseq_pipe.R --expression {params.expression} --species {params.species} \
        --prefix {params.outpre} --signature {params.signature} \
        --outdir {params.outdir} --thread {threads}
        """
if len(ALL_SAMPLES)>1:
    rule scrna_analysis_merge:
        input:
            expression = "Result/QC/%s/%s_filtered_gene_count.h5" % (config["mergedname"],config["mergedname"])
        output:
            specificgene = "Result/Analysis/%s/%s_DiffGenes.tsv" % (config["mergedname"],config["mergedname"]),
            clusterplot = "Result/Analysis/%s/%s_cluster.png" % (config["mergedname"],config["mergedname"]),
            annotateplot = "Result/Analysis/%s/%s_annotated.png" % (config["mergedname"],config["mergedname"]),
            sampleplot = "Result/Analysis/%s/%s_samples.png" %(config["mergedname"],config["mergedname"]),
            tflist = "Result/Analysis/%s/%s.PredictedTFTop10.txt"  % (config["mergedname"],config["mergedname"])
        params:
            expression = "../QC/%s/%s_filtered_gene_count.h5" % (config["mergedname"],config["mergedname"]),
            species = config["species"],
            outpre =  config["mergedname"],
            outdir = "Result/Analysis",
            lisadir = config["lisadir"],
            signature = config["signature"],
            sampleannotation = "../STAR/%s/BarcodeAnnotation.tsv" % config["mergedname"]
        benchmark:
            "Result/Benchmark/%s/%s_Analysis.benchmark" % (config["mergedname"],config["mergedname"])
        threads:
            config["cores"]
        shell:
            """
            python {SCRIPT_PATH}/lisa_path.py --species {params.species} --input {params.lisadir};
            Rscript {RSCRIPT_PATH}/scRNAseq_pipe.R --expression {params.expression} \
            --species {params.species} --prefix {params.outpre} --signature {params.signature} \
            --outdir {params.outdir} --thread {threads} --samplefile {params.sampleannotation}
            """
