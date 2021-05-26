# rule for mapping scATAC sequence

rule scatac_map:
    input:
        fasta = config["genome"]["atac_fasta"],
        r1 = os.path.join("Result/ATAC/Tmp", "%s_R1.barcoded.fastq" %(config["atac_fastqprefix"])),
        r3 = os.path.join("Result/ATAC/Tmp", "%s_R3.barcoded.fastq" %(config["atac_fastqprefix"])),
    output:
        bam = temp("Result/ATAC/minimap2/%s.sortedByPos.bam" %(config["outprefix"])),
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_Minimap2.benchmark" %(config["outprefix"])
    shell:
        "minimap2 -ax sr -t {threads} {input.fasta} {input.r1} {input.r3} "
        "| samtools view --threads {threads} -b "
        "| samtools sort --threads {threads} -o {output.bam};"

rule scatac_fragmentgenerate:
    input:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.bam" %(config["outprefix"]),
    output:
        fragments = "Result/ATAC/minimap2/fragments.tsv",
        bam = "Result/ATAC/minimap2/%s.sortedByPos.CRadded.bam" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/minimap2"
    benchmark:
        "Result/Benchmark/%s_scATAC_FragGenerate.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentGenerate.py -B {input.bam} -O {params.outdir} --addtag CR"

rule scatac_rmdp:
    input:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.CRadded.bam" %(config["outprefix"]),
    output:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.CRadded.rmdp.bam" %(config["outprefix"]),
        metric = "Result/ATAC/minimap2/%s.rmdp.txt" %(config["outprefix"]),
        fragbed = "Result/ATAC/QC/%s_frag.bed" %(config["outprefix"]),
    params:
        sam = "Result/ATAC/minimap2/%s.sortedByPos.CRadded.rmdp.sample.sam" %(config["outprefix"]),
        tmp = "Result/ATAC/Tmp",
    threads:
        int(config["cores"])-2
    benchmark:
        "Result/Benchmark/%s_scATAC_Rmdp.benchmark" %(config["outprefix"])
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR={params.tmp};"
        "samtools view -@ {threads} -s 0.01 -o {params.sam} {input.bam};"
        "awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed};"

rule scatac_bamaddCB:
    input:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.CRadded.rmdp.bam" %(config["outprefix"]),
        bc_correct = "Result/ATAC/minimap2/barcode_correct_uniq.txt"
    output:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/minimap2",
        outprefix = "%s.sortedByPos.rmdp.CBadded" %(config["outprefix"])
    threads:
        int(config["cores"])-2
    benchmark:
        "Result/Benchmark/%s_scATAC_BamAddCB.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_BamAddTag.py -B {input.bam} -T {input.bc_correct} -C CR "
        "-O {params.outdir} -P {params.outprefix};"

