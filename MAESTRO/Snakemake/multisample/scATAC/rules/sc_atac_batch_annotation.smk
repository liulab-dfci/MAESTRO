_annotation_threads = 4

rule scatac_batch_merge_peak_h5:
    input:
        filtercounts = sorted(expand("Result/Analysis/Batch/{sample}/{sample}_filtered_peak_count.h5", sample = ALL_SAMPLES))
    output:
        mergedcount = "Result/Analysis/Batch/all_samples_filtered_peak_count.h5"
    params:
        species = config["species"],
        prefix = sorted(ALL_SAMPLES)
    shell:
        "MAESTRO merge-h5 --type Peak --h5 {input} --species {params.species} --cellprefix {params.prefix} "
        "--directory Result/Analysis/Batch --outprefix all_samples_filtered"


## merge-h5 add the sample_name@ to the cell barcode: sample_name@TGGTCCTTCATTCGGA
rule scatac_batch_merge_gene_score_h5:
    input:
        filter_gene_scores = sorted(expand("Result/Analysis/Batch/{sample}/{sample}_gene_score.h5", sample = ALL_SAMPLES))
    output:
        merged_gene_score = "Result/Analysis/Batch/all_samples_filtered_gene_count.h5"
    params:
        species = config["species"],
        prefix = sorted(ALL_SAMPLES)
    shell:
        "MAESTRO merge-h5 --type Gene --h5 {input} --species {params.species} --cellprefix {params.prefix} "
        "--directory Result/Analysis/Batch --outprefix all_samples_filtered"


rule scatac_batch_merge_fragment:
    input:
        frag = sorted(expand("Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv", sample = ALL_SAMPLES))
    output:
        merged_frag = "Result/Analysis/Batch/all_samples_fragments_corrected_dedup_count.tsv.gz"
    params:
        tmp_fragments = expand("Result/Tmp/{sample}_fragment_tmp.tsv", sample = ALL_SAMPLES)
    shell:
        """
        for i in {input}
            do
                prefix=$(echo $i | sed -E 's,Result/minimap2/(.+)/fragments_corrected_dedup_count.tsv,\\1,')
                echo $prefix
                cat $i | awk -v prefix="$prefix" -v OFS="\t" '$4=prefix"@"$4' > Result/Tmp/${{prefix}}_fragment_tmp.tsv
            done

        cat {params.tmp_fragments} | bgzip > {output}
        rm {params.tmp_fragments}
        tabix -p bed {output}
        """

rule scatac_batch_analysis:
    input:
        filtercount = "Result/Analysis/Batch/all_samples_filtered_peak_count.h5",
        genescore = "Result/Analysis/Batch/all_samples_filtered_gene_count.h5",
        fraggz = "Result/Analysis/Batch/all_samples_fragments_corrected_dedup_count.tsv.gz"
    output:
        specificpeak = "Result/Analysis/Batch/all_samples_DiffPeaks.tsv",
        clusterplot = "Result/Analysis/Batch/all_samples_cluster.png",
        # annotateplot = "Result/Analysis/%s_annotated.png" %(config["outprefix"]),
        tflist = "Result/Analysis/Batch/all_samples.PredictedTFTop10.txt",
        cellcluster = "Result/Analysis/Batch/all_samples_cell_cluster.txt"
    params:
        outdir = "Result/Analysis/Batch",
        genescore = "all_samples_filtered_gene_count.h5",
        outpre = "all_samples",
        count = "all_samples_filtered_peak_count.h5",
        fraggz = "all_samples_fragments_corrected_dedup_count.tsv.gz",
        giggleannotation = config["giggleannotation"],
        species = config["species"],
        signature = config["signature"],
        method = config["method"],
        annotation = config["annotation"],
    threads:
        _annotation_threads
    benchmark:
        "Result/Benchmark/all_samples_Analysis.benchmark" 
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_pipe.R --peakcount {params.count} --rpmatrix {params.genescore} "
        "--species {params.species} --prefix {params.outpre} --annotation {params.annotation} --method {params.method} --signature {params.signature} "
        "--gigglelib {params.giggleannotation} --fragment {params.fraggz} --outdir {params.outdir} --thread {threads}"
