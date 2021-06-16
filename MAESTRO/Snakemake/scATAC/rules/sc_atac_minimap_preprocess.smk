_cat_threads= 2
if config["platform"] == "10x-genomics":
    if config["gzip"]:
        rule scatac_preprocess:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["R1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["R2"],
                r3 = lambda wildcards: FILES[wildcards.sample]["R3"]
            output:
                r1cat = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2cat = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3cat = temp("Result/Tmp/{sample}/{sample}_R3.fastq")
            log:
                "Result/Log/{sample}_preprocess.log"
            benchmark:
                "Result/Benchmark/{sample}_Preprocess.benchmark"
            threads: _cat_threads
            shell:
                """
                gunzip -c {input.r1} > {output.r1cat} 2> {log}
                gunzip -c {input.r2} > {output.r2cat} 2>> {log}
                gunzip -c {input.r3} > {output.r3cat} 2>> {log}
                """
    else:
        rule scatac_link_fastq:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["R1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["R2"],
                r3 = lambda wildcards: FILES[wildcards.sample]["R3"]
            output:
                r1cat = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2cat = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3cat = temp("Result/Tmp/{sample}/{sample}_R3.fastq")
            log:
                "Result/Log/{sample}_preprocess.log"
            benchmark:
                "Result/Benchmark/{sample}_Preprocess.benchmark"
            shell:
                """
                cat {input.r1} > {output.r1cat}
                cat {input.r2} > {output.r2cat}
                cat {input.r3} > {output.r3cat}
                """

elif config["platform"] == "sci-ATAC-seq":
    if config["gzip"]:
        rule scatac_preprocess:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["2"],
            output:
                r1cat = temp("Result/Tmp/{sample}/{sample}_pre_R1.fastq"),
                r2cat = temp("Result/Tmp/{sample}/{sample}_pre_R2.fastq"),
            log:
                "Result/Log/{sample}_preprocess.log"
            benchmark:
                "Result/Benchmark/{sample}_Preprocess.benchmark"
            threads: _cat_threads
            shell:
                """
                gunzip -c {input.r1} > {output.r1cat} 2> {log}
                gunzip -c {input.r2} > {output.r2cat} 2>> {log}
                """

        rule scatac_sci_BarcodeExtract:
            input:
                r1 = "Result/Tmp/{sample}/{sample}_pre_R1.fastq",
                r2 = "Result/Tmp/{sample}/{sample}_pre_R2.fastq",
            output:
                r1 = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2 = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3 = temp("Result/Tmp/{sample}/{sample}_R3.fastq"),
            params:
                outdir = "Result/Tmp/{sample}/"
            benchmark:
                "Result/Benchmark/{sample}_BarcodeExtract.benchmark"
            shell:
                "python " + SCRIPT_PATH + "/scATAC_sci_BarcodeExtract.py --R1 {input.r1} --R2 {input.r2} -O {params.outdir}"

        rule scatac_preprocess:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["2"],
                r3 = lambda wildcards: FILES[wildcards.sample]["3"]
            output:
                r1cat = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2cat = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3cat = temp("Result/Tmp/{sample}/{sample}_R3.fastq")
            log:
                "Result/Log/{sample}_preprocess.log"
            benchmark:
                "Result/Benchmark/{sample}_Preprocess.benchmark"
            threads: _cat_threads
            shell:
                """
                gunzip -c {input.r1} > {output.r1cat} 2> {log}
                gunzip -c {input.r2} > {output.r2cat} 2>> {log}
                gunzip -c {input.r3} > {output.r3cat} 2>> {log}
                """
    else:
        rule scatac_sci_BarcodeExtract:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["2"],
            output:
                r1 = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2 = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3 = temp("Result/Tmp/{sample}/{sample}_R3.fastq"),
            params:
                outdir = "Result/Tmp/{sample}/"
            benchmark:
                "Result/Benchmark/{sample}_BarcodeExtract.benchmark"
            shell:
                "python " + SCRIPT_PATH + "/scATAC_sci_BarcodeExtract.py --R1 {input.r1} --R2 {input.r2} -O {params.outdir}"


        rule scatac_preprocess:
            input:
                r1 = lambda wildcards: FILES[wildcards.sample]["1"],
                r2 = lambda wildcards: FILES[wildcards.sample]["2"],
                r3 = lambda wildcards: FILES[wildcards.sample]["3"],
            output:
                r1cat = temp("Result/Tmp/{sample}/{sample}_R1.fastq"),
                r2cat = temp("Result/Tmp/{sample}/{sample}_R2.fastq"),
                r3cat = temp("Result/Tmp/{sample}/{sample}_R3.fastq"),
            log:
                "Result/Log/{sample}_preprocess.log"
            benchmark:
                "Result/Benchmark/{sample}_Preprocess.benchmark"
            threads: _cat_threads
            shell:
                """
                cat {input.r1} > {output.r1cat}
                cat {input.r2} > {output.r2cat}
                cat {input.r3} > {output.r3cat}
                """

rule scatac_fqaddbarcode:
    input:
        r1 = "Result/Tmp/{sample}/{sample}_R1.fastq",
        r2 = "Result/Tmp/{sample}/{sample}_R2.fastq",
        r3 = "Result/Tmp/{sample}/{sample}_R3.fastq",
    output:
        r1 = temp("Result/Tmp/{sample}/{sample}_R1.barcoded.fastq"),
        r3 = temp("Result/Tmp/{sample}/{sample}_R3.barcoded.fastq"),
    benchmark:
        "Result/Benchmark/{sample}_FqAddbarcode.benchmark"
    shell:
        """
        base=`head -n 2 {input.r2} | tail -n 1 | wc -L`

        sinto barcode --barcode_fastq {input.r2} --read1 {input.r1} --read2 {input.r3} -b $base
        """
