if config["platform"] == "10x-genomics" or config["platform"] == "Dropseq":
    rule scrna_qc:
        input:
            rawmtx = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/matrix.mtx",
            feature = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/features.tsv",
            barcode = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/barcodes.tsv"
        output:
            countgene = "Result/QC/{sample}/{sample}_count_gene_stat.txt",
            filtermatrix = "Result/QC/{sample}/{sample}_filtered_gene_count.h5",
            rnafilterplot = "Result/QC/{sample}/{sample}_scRNA_cell_filtering.png"
        params:
            counts = config["cutoff"]["count"],
            gene = config["cutoff"]["gene"],
            outpre = "{sample}/{sample}",
            outdir = "Result/QC",
            species = config["species"]
        benchmark:
            "Result/Benchmark/{sample}/{sample}_QC.benchmark"
        shell:
            """
            MAESTRO scrna-qc \
                --format mtx \
                --matrix {input.rawmtx} \
                --feature {input.feature} \
                --barcode {input.barcode} \
                --species {params.species} \
                --count-cutoff {params.counts} \
                --gene-cutoff {params.gene} \
                --directory {params.outdir} \
                --outprefix {params.outpre}
            """

    if len(ALL_SAMPLES) > 1:
        rule scrna_qc_merge:
            input:
                rawmtx = "Result/STAR/%s/rawmatrix/matrix.mtx" % config["mergedname"],
                feature = "Result/STAR/%s/rawmatrix/features.tsv" % config["mergedname"],
                barcode = "Result/STAR/%s/rawmatrix/barcodes.tsv" % config["mergedname"]
            output:
                countgene = "Result/QC/%s/%s_count_gene_stat.txt" % (config["mergedname"],config["mergedname"]),
                filtermatrix = "Result/QC/%s/%s_filtered_gene_count.h5" % (config["mergedname"],config["mergedname"]),
                rnafilterplot = "Result/QC/%s/%s_scRNA_cell_filtering.png" % (config["mergedname"],config["mergedname"])
            params:
                counts = config["cutoff"]["count"],
                gene = config["cutoff"]["gene"],
                outpre = config["mergedname"],
                outdir = "Result/QC/%s" % config["mergedname"],
                species = config["species"]
            benchmark:
                "Result/Benchmark/%s/%s_QC.benchmark" % (config["mergedname"],config["mergedname"])
            shell:
                """
                MAESTRO scrna-qc \
                    --format mtx \
                    --matrix {input.rawmtx} \
                    --feature {input.feature} \
                    --barcode {input.barcode} \
                    --species {params.species} \
                    --count-cutoff {params.counts} \
                    --gene-cutoff {params.gene} \
                    --directory {params.outdir} \
                    --outprefix {params.outpre}
                """
elif config["platform"] == "Smartseq2":
    rule scrna_rsem_expre:
        input:
            transbam = "Result/STAR/{sample}/{sample}Aligned.toTranscriptome.out.bam"
        output:
            generesult = "Result/STAR/{sample}/{sample}.genes.results"
        params:
            reference = config["genome"]["rsem"],
            sample = "Result/STAR/{sample}/{sample}"
        log:
            "Result/Log/{sample}/{sample}_RSEM_map.log"
        threads:
            config["cores"]
        shell:
            """
            rsem-calculate-expression \
                --paired-end --bam --estimate-rspd \
                -p {threads} \
                --append-names {input.transbam} {params.reference} {params.sample} \
                > {log} 2>&1
            """

    rule scrna_rsem_count:
        input:
            generesult = "Result/STAR/{sample}/{sample}.genes.results"
        output:
            expression = "Result/Count/{sample}/{sample}_gene_count_matrix.txt"
        shell:
            """
            rsem-generate-data-matrix {input.generesult} > {output.expression}
            """


    rule scrna_qc:
        input:
            expression = "Result/Count/{sample}/{sample}_gene_count_matrix.txt"
        output:
            countgene = "Result/QC/{sample}/{sample}_count_gene_stat.txt",
            filtermatrix = "Result/QC/{sample}/{sample}_filtered_gene_count.h5",
            rnafilterplot = "Result/QC/{sample}/{sample}_scRNA_cell_filtering.png"
        params:
            counts = config["cutoff"]["count"],
            gene = config["cutoff"]["gene"],
            outpre = "{sample}/{sample}",
            outdir = "Result/QC",
            species = config["species"]
        benchmark:
            "Result/Benchmark/{sample}/{sample}_QC.benchmark"
        shell:
            """
            MAESTRO scrna-qc \
                --format plain \
                --matrix {input.expression} \
                --species {params.species} \
                --count-cutoff {params.counts} \
                --gene-cutoff {params.gene} \
                --directory {params.outdir} \
                --outprefix {params.outpre}
            """

    rule scrna_bammerge:
        input:
            genomebam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam"
        output:
            bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            bai = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
            bamlist = "Result/STAR/{sample}/{sample}_bamlist.txt"
        params:
            bamprefix = "Result/STAR/{sample}/{sample}_bamlist_",
            subprefix = "Result/STAR/{sample}/{sample}"
        benchmark:
            "Result/Benchmark/{sample}/{sample}_BamMerge.benchmark"
        threads:
            config["cores"]
        shell:
            """
            ls Result/STAR/*Aligned.sortedByCoord.out.bam > {output.bamlist};
            split -1000 -d {output.bamlist} {params.bamprefix};
            for file in $(ls {params.bamprefix}*);
            do
                sub=${{file#{params.bamprefix}}};
                samtools merge \
                    --threads {threads} \
                    {params.subprefix}.${{sub}}.Aligned.sortedByCoord.out.bam \
                    -b ${{file}};
            done;
            samtools merge \
                --threads {threads} \
                {output.bam} \
                {params.subprefix}.*.Aligned.sortedByCoord.out.bam;
            rm {params.subprefix}.[0-9]*.Aligned.sortedByCoord.out.bam;
            samtools index -b \
                -@ {threads} \
                {output.bam}
            """
