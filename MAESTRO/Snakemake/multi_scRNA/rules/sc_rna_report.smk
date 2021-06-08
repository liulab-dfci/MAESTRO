
if config["rseqc"]:
    rule scrna_report:
        input:
            clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
            annotateplot = "Result/Analysis/{sample}/{sample}_annotated.png",
            tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt",
            readdistrplot = "Result/QC/{sample}/{sample}_scRNA_read_distr.png",
            qualplot = "Result/QC/{sample}/{sample}_scRNA_read_quality.png",
            nvcplot = "Result/QC/{sample}/{sample}_scRNA_NVC.png",
            gcplot = "Result/QC/{sample}/{sample}_scRNA_GCcontent.png",
            genecovplot = "Result/QC/{sample}/{sample}_scRNA_genebody_cov.png",
            rnafilterplot = "Result/QC/{sample}/{sample}_scRNA_cell_filtering.png",
        output:
            summaryreport = "Result/{sample}/{sample}_scRNA_report.html",
        params:
            outdir = "Result",
            outpre = "{sample}/{sample}",
            fastqdir = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
            species = config["species"],
            platform = config["platform"],
        benchmark:
            "Result/Benchmark/{sample}/{sample}_Report.benchmark"
        shell:
            """
			python {SCRIPT_PATH}/scRNA_HTMLReport.py \
			--directory {params.outdir} \
			--outprefix {params.outpre} \
            --fastq-dir {params.fastqdir} \
			--species {params.species} \
			--platform {params.platform} \
			--rseqc
			"""
else:
    rule scrna_report:
        input:
            clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
            annotateplot = "Result/Analysis/{sample}/{sample}_annotated.png",
            tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt",
            rnafilterplot = "Result/QC/{sample}/{sample}_scRNA_cell_filtering.png",
        output:
            summaryreport = "Result/{sample}/{sample}_scRNA_report.html",
        params:
            outdir = "Result",
            outpre = "{sample}/{sample}",
            fastqdir = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
            species = config["species"],
            platform = config["platform"],
            rseqc = False
        benchmark:
            "Result/Benchmark/{sample}/{sample}_Report.benchmark"
        shell:
            """
			python {SCRIPT_PATH}/scRNA_HTMLReport.py \
			--directory {params.outdir} \
			--outprefix {params.outpre} \
            --fastq-dir {params.fastqdir} \
			--species {params.species} \
			--platform {params.platform}
			"""

if len(ALL_SAMPLES) > 1:
    rule scrna_report_merge:
        input:
            clusterplot = "Result/Analysis/%s_cluster.png" % config["mergedname"],
            sampleplot = "Result/Analysis/%s_samples.png" %config["mergedname"],
            annotateplot = "Result/Analysis/%s_annotated.png" % config["mergedname"],
            tflist = "Result/Analysis/%s.PredictedTFTop10.txt" % config["mergedname"],
            rnafilterplot = "Result/QC/%s_scRNA_cell_filtering.png" % config["mergedname"]
        output:
            summaryreport = "Result/%s_scRNA_report.html" % config["mergedname"]
        params:
            outdir = "Result",
            outpre = config["mergedname"],
            species = config["species"],
            platform = config["platform"],
            rseqc = False
        benchmark:
            "Result/Benchmark/%s_Report.benchmark" % config["mergedname"]
        shell:
            """
			python {SCRIPT_PATH}/scRNA_HTMLReport.py \
			--directory {params.outdir} \
			--outprefix {params.outpre} \
			--species {params.species} \
			--platform {params.platform} \
            --multisample
			"""