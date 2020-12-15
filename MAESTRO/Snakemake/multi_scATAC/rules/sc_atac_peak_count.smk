_countpeak_threads = 4

rule scatac_countpeak:
    input:
        finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed",
        validbarcode = "Result/QC/{sample}/{sample}_scATAC_validcells.txt",
        frag = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv"
    output:
        counts = "Result/Analysis/{sample}/{sample}_peak_count.h5"
    params:
        species = config["species"],
        outdir = "Result/Analysis/{sample}",
        outpre = "{sample}"
    threads:
        _countpeak_threads
    benchmark:
        "Result/Benchmark/{sample}_PeakCount.benchmark"
    shell:
        """
        MAESTRO scatac-peakcount --peak {input.finalpeak} --fragment {input.frag} --barcode {input.validbarcode} \
        --species {params.species} --cores {threads} --directory {params.outdir} --outprefix {params.outpre}
        """

rule scatac_qcfilter:
    input:
        counts = "Result/Analysis/{sample}/{sample}_peak_count.h5"
    output:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5"
    params:
        outdir = "Result/QC/{sample}",
        outpre = "{sample}",
        peak = config["cutoff"]["peak"],
        cell = config["cutoff"]["cell"]
    benchmark:
        "Result/Benchmark/{sample}_QCFilter.benchmark"
    shell:
        """
        MAESTRO scatac-qc --format h5 --peakcount {input.counts} --peak-cutoff {params.peak} --cell-cutoff {params.cell} \
        --directory {params.outdir} --outprefix {params.outpre}
        """

rule scatac_genescore:
    input:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5",
        genebed = "%s/annotations/%s_ensembl.bed" %(SCRIPT_PATH, config["species"]),
    output:
        genescore = "Result/Analysis/{sample}/{sample}_gene_score.h5"
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/Analysis/{sample}",
        outpre = "{sample}",
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/{sample}_GeneScore.benchmark"
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"
