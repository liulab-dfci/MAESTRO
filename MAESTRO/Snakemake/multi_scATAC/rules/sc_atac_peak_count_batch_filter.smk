rule scatac_batch_qcfilter:
    input:
        counts = "Result/Analysis/Batch/{sample}/{sample}_peak_count.h5"
    output:
        filtercount = "Result/Analysis/Batch/{sample}/{sample}_filtered_peak_count.h5"
    params:
        outdir = "Result/Analysis/Batch/{sample}",
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

rule scatac_batch_genescore:
    input:
        filtercount = "Result/Analysis/Batch/{sample}/{sample}_filtered_peak_count.h5",
        genebed = "%s/annotations/%s_ensembl.bed" %(SCRIPT_PATH, config["species"]),
    output:
        genescore = "Result/Analysis/Batch/{sample}/{sample}_gene_score.h5"
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/Analysis/Batch/{sample}",
        outpre = "{sample}",
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/{sample}_GeneScore.benchmark"
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"
