
_annotation_threads = 4

## this is too verbose. scATACseq_pipe.R setwd() which is not good. see we can improve it.
def get_fragments_zip1(wildcards):
    if config["deduplication"] == "cell level":
        return "../../minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv.gz".format(sample = wildcards.sample)
    elif config["deduplication"] == "bulk level":
        return "../../minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv.gz".format(sample = wildcards.sample)
    else:
        print("please specify 'cell level' or 'bulk level")
        sys.exit(1)

def get_fragments_zip2(wildcards):
    if config["deduplication"] == "cell level":
        return "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv.gz".format(sample = wildcards.sample)
    elif config["deduplication"] == "bulk level":
        return "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv.gz".format(sample = wildcards.sample)
    else:
        print("please specify 'cell level' or 'bulk level")
        sys.exit(1)

rule scatac_analysis:
    input:
        filtercount = "Result/QC/{sample}/{sample}_filtered_peak_count.h5",
        genescore = "Result/Analysis/{sample}/{sample}_gene_score.h5",
        fraggz = get_fragments_zip2
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
        count = "../../QC/{sample}/{sample}_filtered_peak_count.h5",
        fraggz = get_fragments_zip1,
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
        "Rscript " + RSCRIPT_PATH + "/scATACseq_pipe.R --peakcount {params.count} --rpmatrix {params.genescore} "
        "--species {params.species} --prefix {params.outpre} --annotation {params.annotation} --method {params.method} --signature {params.signature} "
        "--gigglelib {params.giggleannotation} --fragment {params.fraggz} --outdir {params.outdir} --thread {threads}"
