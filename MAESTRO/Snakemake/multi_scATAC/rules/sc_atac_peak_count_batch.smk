
## when one processes multiple samples from the same experiment, he/she
## will want to merge all the peaks across samples together and then get a new count matrix
## with the new peak set.

_countpeak_threads = 4

## consensus peak calling.
rule scatac_merge_peaks_batch:
    input: expand("Result/Analysis/{sample}/{sample}_final_peaks.bed", sample = ALL_SAMPLES)
    output: "Result/Analysis/Batch/all_samples_peaks.bed"
    params:
        catpeaksort = "Result/Analysis/Batch/all_samples_peaks_cat.bed",
        catpeakbdg = "Result/Analysis/Batch/all_samples_peaks_cat.bdg",
        cutoff_samples = config["cutoff_samples"],
        chrom_len = config["chrom_len"]
    log: "Result/Log/merge_peaks_batch.log"
    shell:
        """
        # based on Tao Liu's gist https://gist.github.com/taoliu/5772f4b44767ff07b85ba11c42a57b78
        # define the minlen and maxgap for peak calling (arbitrary)
        MINLEN=200
        MAXGAP=30

        for f in {input}; do
            bedtools sort -i $f | mergeBed -i - | cut -f 1-3 | grep -v '_' | grep -v 'chrEBV' >> {params.catpeaksort}
        done

        bedtools sort -i {params.catpeaksort} | bedtools genomecov -bga -i - -g {params.chrom_len} > {params.catpeakbdg}

        macs2 bdgpeakcall -i {params.catpeakbdg} -o {output} --no-trackline -c {params.cutoff_samples} -g $MAXGAP -l $MINLEN
        rm {params.catpeaksort}
        rm {params.catpeakbdg}

        """

rule scatac_countpeak_batch:
    input:
        finalpeak = "Result/Analysis/Batch/all_samples_peaks.bed",
        validbarcode = "Result/QC/{sample}/{sample}_scATAC_validcells.txt",
        frag = get_fragments
    output:
        counts = "Result/Analysis/Batch/{sample}/{sample}_peak_count.h5"
    params:
        species = config["species"],
        outdir = "Result/Analysis/Batch/{sample}",
        outpre = "{sample}"
    threads:
        _countpeak_threads
    benchmark:
        "Result/Benchmark/{sample}_PeakCount_batch.benchmark"
    shell:
        """
        MAESTRO scatac-peakcount --peak {input.finalpeak} --fragment {input.frag} --barcode {input.validbarcode} \
        --species {params.species} --cores {threads} --directory {params.outdir} --outprefix {params.outpre}
        """
