_picard_threads = 4
rule scatac_fragmentgenerate:
        input:
            bam = lambda wildcards: FILES[wildcards.sample]
        output:
            fragments = "Result/Mapping/{sample}/fragments.tsv",
        params:
            outdir = "Result/Mapping/{sample}"
        benchmark:
            "Result/Benchmark/{sample}_FragGenerate.benchmark"
        shell:
            "python " + SCRIPT_PATH + "/scATAC_FragmentGenerate.py -B {input.bam} -O {params.outdir} --CBtag CB --count;"

rule scatac_rmdp:
    input:
        bam = lambda wildcards: FILES[wildcards.sample],
    output:
        bam = "Result/Mapping/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam",
        metric = "Result/Mapping/{sample}/{sample}.rmdp.txt",
        fragbed = "Result/QC/{sample}/{sample}_frag.bed"
    params:
        sam = "Result/Mapping/{sample}/{sample}.sortedByPos.CRadded.rmdp.sample.sam"
    threads:
        _picard_threads
    benchmark:
        "Result/Benchmark/{sample}_Rmdp.benchmark"
    shell:
        """
        picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR=Result/Tmp;
        samtools view -@ {threads} -s 0.01 -o {params.sam} {input.bam};
        awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed}
        """

rule scatac_bamindex:
    input:
        bam = "Result/Mapping/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam",
    output:
        bai = "Result/Mapping/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam.bai",
    threads:
        _bamindex_threads
    benchmark:
        "Result/Benchmark/{sample}_BamIndex.benchmark"
    shell:
        "samtools index -@ {threads} {input.bam}"
