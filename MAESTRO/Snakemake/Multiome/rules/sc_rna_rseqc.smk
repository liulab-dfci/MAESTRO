# rule for bulk-level quality control (rseqc)

rule scrna_samplebam:
    input:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
    output:
        bamsample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam" %(config["outprefix"]),
        baisample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam.bai" %(config["outprefix"]),
    benchmark:
        "Result/Benchmark/%s_scRNA_BamSample.benchmark" %(config["outprefix"])
    threads:
        config["cores"]
    shell:
        "samtools view -@ {threads} -s 0.01 -b -o {output.bamsample} {input.bam};"
        "samtools index -@ {threads} {output.bamsample}"

rule scrna_rseqc_readqual:
    input:
        bamsample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam" %(config["outprefix"]),
    output:
        qual = "Result/RNA/QC/%s.qual.r" %(config["outprefix"])
    params:
        outdirpre = "Result/RNA/QC/%s" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcReadqual.benchmark" %(config["outprefix"])
    shell:
        "read_quality.py -i {input.bamsample} -o {params.outdirpre}"

rule scrna_rseqc_nvc:
    input:
        bamsample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam" %(config["outprefix"]),
    output:
        nvc = "Result/RNA/QC/%s.NVC.xls" %(config["outprefix"])
    params:
        outdirpre = "Result/RNA/QC/%s" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcNVC.benchmark" %(config["outprefix"])
    shell:
        "read_NVC.py -i {input.bamsample} -o {params.outdirpre}"

rule scrna_rseqc_gc:
    input:
        bamsample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam" %(config["outprefix"]),
        baisample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam.bai" %(config["outprefix"]),
    output:
        gc = "Result/RNA/QC/%s.GC.xls" %(config["outprefix"])
    params:
        outdirpre = "Result/RNA/QC/%s" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcGC.benchmark" %(config["outprefix"])
    shell:
        "read_GC.py -i {input.bamsample} -o {params.outdirpre}"

rule scrna_rseqc_bamstat:
    input:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
        bai = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam.bai" %(config["outprefix"]),
    output:
        stat = "Result/RNA/QC/%s_bam_stat.txt" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcBamstat.benchmark" %(config["outprefix"])
    shell:
        "bam_stat.py -i {input.bam} > {output.stat};"

rule scrna_rseqc_distr:
    input:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
        bai = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam.bai" %(config["outprefix"]),
        genome = "%s/annotations/%s_RefSeq.bed" %(SCRIPT_PATH, config["species"]),
    output:
        distr = "Result/RNA/QC/%s_read_distribution.txt" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcDistr.benchmark" %(config["outprefix"])
    shell:
        "read_distribution.py -i {input.bam} -r {input.genome} > {output.distr}"

rule scrna_rseqc_genecov:
    input:
        bamsample = "Result/RNA/STAR/%sAligned.sortedByCoord.out.sample.bam" %(config["outprefix"]),
        hkgene = "%s/annotations/%s_HouseKeepingGenes.bed" %(SCRIPT_PATH, config["species"]),
    output:
        genecov = "Result/RNA/QC/%s.geneBodyCoverage.txt" %(config["outprefix"]),
    params:
        outdirpre = "Result/RNA/QC/%s" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_RseqcGenecov.benchmark" %(config["outprefix"])
    shell:
        "geneBody_coverage.py -r {input.hkgene} -i {input.bamsample} -o {params.outdirpre}"

rule scrna_rseqc_plot:
    input:
        stat = "Result/RNA/QC/%s_bam_stat.txt" %(config["outprefix"]),
        distr = "Result/RNA/QC/%s_read_distribution.txt" %(config["outprefix"]),
        qual = "Result/RNA/QC/%s.qual.r" %(config["outprefix"]),
        nvc = "Result/RNA/QC/%s.NVC.xls" %(config["outprefix"]),
        gc = "Result/RNA/QC/%s.GC.xls" %(config["outprefix"]),
        genecov = "Result/RNA/QC/%s.geneBodyCoverage.txt" %(config["outprefix"]),
        # countgene = "Result/QC/%s_count_gene_stat.txt" %(config["outprefix"]),
    output:
        readdistrplot = "Result/RNA/QC/%s_scRNA_read_distr.png" %(config["outprefix"]),
        qualplot = "Result/RNA/QC/%s_scRNA_read_quality.png" %(config["outprefix"]),
        nvcplot = "Result/RNA/QC/%s_scRNA_NVC.png" %(config["outprefix"]),
        gcplot = "Result/RNA/QC/%s_scRNA_GCcontent.png" %(config["outprefix"]),
        genecovplot = "Result/RNA/QC/%s_scRNA_genebody_cov.png" %(config["outprefix"]),
    params:
        outpre = config["outprefix"],
        outdir = "Result/RNA/QC",
        rseqc = "TRUE",
        stat = "%s_bam_stat.txt" %(config["outprefix"]),
        distr = "%s_read_distribution.txt" %(config["outprefix"]),
        qual = "%s.qual.r" %(config["outprefix"]),
        nvc = "%s.NVC.xls" %(config["outprefix"]),
        gc = "%s.GC.xls" %(config["outprefix"]),
        genecov = "%s.geneBodyCoverage.txt" %(config["outprefix"]),
        # countgene = "%s_count_gene_stat.txt" %(config["outprefix"]),
        # count = config["cutoff"]["count"],
        # gene = config["cutoff"]["gene"]
    benchmark:
        "Result/Benchmark/%s_scRNA_QCPlot.benchmark" %(config["outprefix"])
    shell:
        "Rscript " + RSCRIPT_PATH + "/scRNAseq_qc.R --prefix {params.outpre} --outdir {params.outdir} "
        "--bamstat {params.stat} --readdistr {params.distr} --qual {params.qual} --nvc {params.nvc} "
        "--gc {params.gc} --genecov {params.genecov}"
