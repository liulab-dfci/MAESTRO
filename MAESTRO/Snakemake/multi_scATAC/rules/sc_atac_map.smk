_minmap_threads = 8
_picard_threads = 4
_bamAddCB_threads = 4
_bamindex_threads = 4

rule scatac_map:
    input:
        fasta = config["genome"]["fasta"],
        r1 = "Result/Tmp/{sample}/{sample}_R1.barcoded.fastq",
        r3 = "Result/Tmp/{sample}/{sample}_R3.barcoded.fastq"
    output:
            bam = temp("Result/minimap2/{sample}/{sample}.sortedByPos.bam")
    threads:
        _minmap_threads
    benchmark:
            "Result/Benchmark/{sample}_Minimap2.benchmark"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.r1} {input.r3} \
        | samtools view --threads {threads} -b \
        | samtools sort --threads {threads} -o {output.bam}
        """

rule scatac_fragmentgenerate:
        input:
            bam = "Result/minimap2/{sample}/{sample}.sortedByPos.bam"
        output:
            fragments = "Result/minimap2/{sample}/fragments.tsv",
            bam = "Result/minimap2/{sample}/{sample}.sortedByPos.CRadded.bam"
        params:
            outdir = "Result/minimap2/{sample}"
        benchmark:
            "Result/Benchmark/{sample}_FragGenerate.benchmark" 
        shell:
            "python " + SCRIPT_PATH + "/scATAC_FragmentGenerate.py -B {input.bam} -O {params.outdir} --addtag CR"

rule scatac_rmdp:
    input:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.CRadded.bam",
    output:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.CRadded.rmdp.bam",
        metric = "Result/minimap2/{sample}/{sample}.rmdp.txt",
        fragbed = "Result/QC/{sample}/{sample}_frag.bed"
    params:
        sam = "Result/minimap2/{sample}/{sample}.sortedByPos.CRadded.rmdp.sample.sam"
    threads:
        _picard_threads
    benchmark:
        "Result/Benchmark/{sample}_Rmdp.benchmark" 
    shell:
        """
        picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR=Result/Tmp

        samtools view -@ {threads} -s 0.01 -o {params.sam} {input.bam}

        awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed}
        """

rule scatac_bamaddCB:
        input:
            bam = "Result/minimap2/{sample}/{sample}.sortedByPos.CRadded.rmdp.bam",
            bc_correct = "Result/minimap2/{sample}/barcode_correct_uniq.txt"
        output:
            bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam"
        params:
            outdir = "Result/minimap2/{sample}",
            outprefix = "{sample}.sortedByPos.rmdp.CBadded" 
        threads:
            _bamAddCB_threads
        benchmark:
            "Result/Benchmark/{sample}_BamAddCB.benchmark"
        shell:
            "python " +  SCRIPT_PATH + "/scATAC_BamAddTag.py -B {input.bam} -T {input.bc_correct} -C CR "
            "-O {params.outdir} -P {params.outprefix}"

            


rule scatac_bamindex:
    input:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam",
    output:
        bai = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam.bai",
    threads:
        _bamindex_threads
    benchmark:
        "Result/Benchmark/{sample}_BamIndex.benchmark" 
    shell:
        "samtools index -@ {threads} {input.bam}"




