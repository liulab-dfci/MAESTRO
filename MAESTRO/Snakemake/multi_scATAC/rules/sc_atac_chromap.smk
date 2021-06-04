_chromap_threads = 8

rule scatac_preprocess:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]["R1"],
        r2 = lambda wildcards: FILES[wildcards.sample]["R2"],
        r3 = lambda wildcards: FILES[wildcards.sample]["R3"]
    output:
        r1cat = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
        r2cat = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
        r3cat = temp("Result/Tmp/{sample}/{sample}_R3.fastq"),
    log:
        "Result/Log/{sample}_preprocess.log"
    benchmark:
        "Result/Benchmark/{sample}_Preprocess.benchmark"
    shell:
        """
        gunzip -c {input.r1} > {output.r1cat} 2> {log}
        gunzip -c {input.r2} > {output.r2cat} 2>> {log}
        gunzip -c {input.r3} > {output.r3cat} 2>> {log}
        """

rule scatac_chromap:
    input:
        r1 = "Result/Tmp/{sample}/{sample}_R1.fastq",
        r3 = "Result/Tmp/{sample}/{sample}_R3.fastq",
        barcode = "Result/Tmp/{sample}/{sample}_R2.fastq",
    output:
        "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv",
    benchmark:
        "Result/Benchmark/{sample}_chromap.benchmark"
    threads:
        _chromap_threads
    params:
        index = config["genome"]["index"],
        fasta = config["genome"]["fasta"],
        whitelist = config["whitelist"],
    shell:
        "chromap --preset atac -x {params.index} -r {params.fasta} -1 {input.r1} -2 {input.r3} -o {output} -b {input.barcode} -t {threads} --barcode-whitelist {params.whitelist}"


rule scatac_process_frag:
    input:
        "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv"
    output:
        fraggz = "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv.gz"
    shell:
        "bgzip -c {input} > {output.fraggz};"
        "tabix -p bed {output.fraggz}"
