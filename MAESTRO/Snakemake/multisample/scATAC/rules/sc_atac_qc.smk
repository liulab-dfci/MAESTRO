_samtools_thead = 4
_qcplot_threads = 4

rule scatac_qcstat_mapped:
    input:
        frag_count = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv",
    output:
        mapped_stat = temp("Result/QC/{sample}/singlecell_mapped.txt"),
    params:
        frag_count_sort = "Result/minimap2/{sample}/fragments_corrected_count_sortedbybarcode.tsv"
    benchmark:
        "Result/Benchmark/{sample}_SingleQCMappability.benchmark" 
    shell:
        """
        sort -k4,4 -V {input.frag_count} > {params.frag_count_sort}

        bedtools groupby -i {params.frag_count_sort} -g 4 -c 5 -o sum > {output.mapped_stat}
        """


rule scatac_qcstat_promoter:
    input:
        frag_count = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv",
        promoter = config["promoter"]
    output:
        promoter_stat = temp("Result/QC/{sample}/singlecell_promoter.txt"),
    params:
        frag_promoter = "Result/minimap2/{sample}/fragments_promoter.tsv",
        frag_promoter_sort = "Result/minimap2/{sample}/fragments_promoter_sortbybarcode.tsv",
    benchmark:
        "Result/Benchmark/{sample}_SingleQCPromoter.benchmark" 
    shell:
        """
        bedtools intersect -wa -a {input.frag_count} -b {input.promoter} -u > {params.frag_promoter}

        sort -k4,4 -V {params.frag_promoter} > {params.frag_promoter_sort}

        bedtools groupby -i {params.frag_promoter_sort} -g 4 -c 5 -o sum > {output.promoter_stat}
        """


rule scatac_qcstat_singlecell:
    input:
        mapped_stat = "Result/QC/{sample}/singlecell_mapped.txt",
        promoter_stat = "Result/QC/{sample}/singlecell_promoter.txt",
    output:
        single_stat = "Result/QC/{sample}/singlecell.txt"
    params:
        mapped_stat_sort = "Result/QC/{sample}/singlecell_mapped_sortbybarcode.txt",
        promoter_stat_sort = "Result/QC/{sample}/singlecell_promoter_sortbybarcode.txt"
    benchmark:
        "Result/Benchmark/{sample}_SingleQCStat.benchmark" 
    shell:
        """
        sort -k1,1 {input.mapped_stat} > {params.mapped_stat_sort}

        sort -k1,1 {input.promoter_stat} > {params.promoter_stat_sort}

        join --nocheck-order -t $'\t' -a1 -e'0' -o'1.1 1.2 2.2' -1 1 -2 1 {params.mapped_stat_sort} {params.promoter_stat_sort} > {output.single_stat}
        """



rule scatac_qcstat_bulk:
    input:
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.bam",
        promoter = config["promoter"],
        peak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak"
    output:
        bulk_stat = "Result/QC/{sample}/flagstat.txt",
        bam = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.unique.bam",
        bed = "Result/minimap2/{sample}/{sample}.sortedByPos.rmdp.CBadded.unique.bed",
    threads:
        _samtools_thead
    benchmark:
        "Result/Benchmark/{sample}/{sample}_BulkQCStat.benchmark"
    shell:
        """
        samtools flagstat --threads {threads} {input.bam} > {output.bulk_stat}
        samtools view -F 2316 -f 0x2 -q 30 -b -o {output.bam} {input.bam}
        samtools view {output.bam} -c >> {output.bulk_stat}
        bedtools bamtobed -i {output.bam} > {output.bed}
        grep 'chrM' {output.bed} -c >> {output.bulk_stat} || true
        grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.bulk_stat} || true
        grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.bulk_stat} || true 
        """



rule scatac_qcplot:
    input:
        fragbed = "Result/QC/{sample}/{sample}_frag.bed",
        single_stat = "Result/QC/{sample}/singlecell.txt",
        bulk_stat = "Result/QC/{sample}/flagstat.txt"
    output:
        readdistr = "Result/QC/{sample}/{sample}_scATAC_read_distr.png",
        qcfrag = "Result/QC/{sample}/{sample}_scATAC_fragment_size.png",
        qcfrip = "Result/QC/{sample}/{sample}_scATAC_cell_filtering.png",
        validbarcode = "Result/QC/{sample}/{sample}_scATAC_validcells.txt"
    params:
        outdir = "Result/QC/{sample}",
        outpre = "{sample}",
        fragbed = "{sample}_frag.bed",
        single_stat = "singlecell.txt",
        bulk_stat = "flagstat.txt",
        count = config["cutoff"]["count"],
        frip = config["cutoff"]["frip"]
    threads:
        _qcplot_threads
    benchmark:
        "Result/Benchmark/{sample}_QCPlot.benchmark" 
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_qc.R --bulkstat {params.bulk_stat} --fragment {params.fragbed} --singlestat {params.single_stat} "
        "--countcutoff {params.count} --fripcutoff {params.frip} --prefix {params.outpre} --outdir {params.outdir}"
        

