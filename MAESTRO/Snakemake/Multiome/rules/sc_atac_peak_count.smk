# rule for counting peaks

rule scatac_countpeak:
    input:
        finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"]),
        validbarcode = "Result/ATAC/QC/%s_scATAC_validcells.txt" %(config["outprefix"]),
        frag = "Result/ATAC/minimap2/fragments_corrected_count.tsv"
    output:
        counts = "Result/ATAC/Analysis/" + config["outprefix"] + "_peak_count.h5"
    params:
        species = config["species"],
        outdir = "Result/ATAC/Analysis",
        outpre = config["outprefix"]
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_PeakCount.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO scatac-peakcount --binary --peak {input.finalpeak} --fragment {input.frag} --barcode {input.validbarcode} "
        "--species {params.species} --cores {threads} --directory {params.outdir} --outprefix {params.outpre}"

rule scatac_qcfilter:
    input:
        counts = "Result/ATAC/Analysis/%s_peak_count.h5" %(config["outprefix"]),
    output:
        filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/QC",
        outpre = config["outprefix"],
        peak = config["cutoff"]["peak"],
        cell = config["cutoff"]["cell"],
    benchmark:
        "Result/Benchmark/%s_scATAC_QCFilter.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO scatac-qc --format h5 --peakcount {input.counts} --peak-cutoff {params.peak} --cell-cutoff {params.cell} "
        "--directory {params.outdir} --outprefix {params.outpre}"

rule scatac_genescore:
    input:
        filtercount = "Result/ATAC/QC/%s_filtered_peak_count.h5" %(config["outprefix"]),
        genebed = "%s/annotations/%s_ensembl.bed" %(SCRIPT_PATH, config["species"])
    output:
        genescore = "Result/ATAC/Analysis/%s_gene_score.h5" %(config["outprefix"])
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/ATAC/Analysis",
        outpre = config["outprefix"],
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/%s_scATAC_GeneScore.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"

