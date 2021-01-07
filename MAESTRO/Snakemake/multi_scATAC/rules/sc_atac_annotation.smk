
_annotation_threads = 4

## this is too verbose. scATACseq_pipe.R setwd() which is not good. see we can improve it.

rule scatac_analysis:
    input:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5",
        genescore = "Result/Analysis/{sample}/{sample}_gene_score.h5",
        fraggz = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv.gz"
    output:
        specificpeak = "Result/Analysis/{sample}/{sample}_DiffPeaks.tsv",
        clusterplot = "Result/Analysis/{sample}/{sample}_cluster.png",
        # annotateplot = "Result/Analysis/%s_annotated.png" %(config["outprefix"]),
        tflist = "Result/Analysis/{sample}/{sample}.PredictedTFTop10.txt",
        cellcluster = "Result/Analysis/{sample}/{sample}_cell_cluster.txt"
    params:
        outdir = "Result/Analysis/{sample}",
        genescore = "{sample}_gene_score.h5",
        outpre = "{sample}",
        counts = "../../QC/{sample}/{sample}_filtered_peak_count.h5",
        fraggz = "../../minimap2/{sample}/fragments_corrected_dedup_count.tsv.gz",
        giggleannotation = config["giggleannotation"],
        species = config["species"],
        signature = config["signature"],
        method = config["method"],
        annotation = config["annotation"],
    threads:
        _annotation_threads
    benchmark:
        "Result/Benchmark/{sample}_Analysis.benchmark"
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_pipe.R --peakcount {params.counts} --rpmatrix {params.genescore} "
        "--species {params.species} --prefix {params.outpre} --annotation {params.annotation} --method {params.method} --signature {params.signature} "
        "--gigglelib {params.giggleannotation} --fragment {params.fraggz} --outdir {params.outdir} --thread {threads}"
