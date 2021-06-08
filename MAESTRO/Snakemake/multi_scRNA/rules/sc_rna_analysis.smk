
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
            expression = "Result/QC/%s_filtered_gene_count.h5" % config["mergedname"]
        output:
            specificgene = "Result/Analysis/%s_DiffGenes.tsv" % config["mergedname"],
            clusterplot = "Result/Analysis/%s_cluster.png" % config["mergedname"],
            annotateplot = "Result/Analysis/%s_annotated.png" % config["mergedname"],
            sampleplot = "Result/Analysis/%s_samples.png" %config["mergedname"],
            tflist = "Result/Analysis/%s.PredictedTFTop10.txt"  % config["mergedname"]   
        params:
            expression = "../QC/%s_filtered_gene_count.h5" % config["mergedname"],
            species = config["species"],
            outpre =  config["mergedname"],
            outdir = "Result/Analysis",
            lisadir = config["lisadir"],
            signature = config["signature"],
            sampleannotation = "../%s/BarcodeAnnotation.tsv" % config["mergedname"]
        benchmark:
            "Result/Benchmark/%s_Analysis.benchmark" % config["mergedname"]
        threads:
            config["cores"]
        shell:
            """
            python {SCRIPT_PATH}/lisa_path.py --species {params.species} --input {params.lisadir};
            Rscript {RSCRIPT_PATH}/scRNAseq_pipe.R --expression {params.expression} \
            --species {params.species} --prefix {params.outpre} --signature {params.signature} \
            --outdir {params.outdir} --thread {threads} --samplefile {params.sampleannotation}
            """