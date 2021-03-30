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


if config["rseqc"]:
    rule multiome_report:
        input:
            mergeobject = "Result/Multiome/%s_multiome_Object.rds" %(config["outprefix"]),
            readdistrplot = "Result/RNA/QC/%s_scRNA_read_distr.png" %(config["outprefix"]),
        output:
            summaryreport = "Result/%s_multiome_report.html" %(config["outprefix"]),
        params:
            outdir = "Result",
            outpre = config["outprefix"],
            atac_fastqdir = config["atac_fastqdir"],
            rna_fastqdir = config["rna_fastqdir"],
            species = config["species"],
        benchmark:
            "Result/Benchmark/%s_Report.benchmark" %(config["outprefix"])
        shell:
            "python " + SCRIPT_PATH + "/Multiome_HTMLReport.py --directory {params.outdir} --outprefix {params.outpre} "
            "--rna-fastq-dir {params.rna_fastqdir} --atac-fastq-dir {params.atac_fastqdir} --species {params.species}"

else:
    rule multiome_report:
        input:
            mergeobject = "Result/Multiome/%s_multiome_Object.rds" %(config["outprefix"]),
        output:
            summaryreport = "Result/%s_multiome_report.html" %(config["outprefix"]),
        params:
            outdir = "Result",
            outpre = config["outprefix"],
            atac_fastqdir = config["atac_fastqdir"],
            rna_fastqdir = config["rna_fastqdir"],
            species = config["species"],
        benchmark:
            "Result/Benchmark/%s_Report.benchmark" %(config["outprefix"])
        shell:
            "python " + SCRIPT_PATH + "/Multiome_HTMLReport.py --directory {params.outdir} --outprefix {params.outpre} "
            "--rna-fastq-dir {params.rna_fastqdir} --atac-fastq-dir {params.atac_fastqdir} --species {params.species}"
