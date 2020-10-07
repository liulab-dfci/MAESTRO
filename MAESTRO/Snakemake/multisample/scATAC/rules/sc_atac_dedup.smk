
## this function is used for downstream workflow (peak call, peak count, qc) to determine if the input fragment file is deduplicated
## at cell level or bulk level

import sys

def get_fragments(wildcards):
    if config["deduplication"] == "cell level":
        return "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv".format(sample = wildcards.sample)
    elif config["deduplication"] == "bulk level":
        return "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv".format(sample = wildcards.sample)
    else:
        print("please specify 'cell level' or 'bulk level")
        sys.exit(1)

if config["deduplication"] == "cell level":
    rule scatac_celldedup:
        input:
            frag_count = "Result/minimap2/{sample}/fragments_corrected_count.tsv",
        output:
            frag_dedup = "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv",
        params:
            source_dir = os.path.dirname(os.path.dirname(srcdir("Snakefile"))) # strip the /rules from srcdir("Snakefile") 
        benchmark:
            "Result/Benchmark/{sample}_Celldedup.benchmark" 
        shell:
            "ln -s {params.source_dir}/{input.frag_count} {output.frag_dedup};"
elif config["deduplication"] == "bulk level":
    rule scatac_bulkdedup:
        input:
            frag_count = "Result/minimap2/{sample}/fragments_corrected_count.tsv",
        output:
            frag_dedup = "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv",
        benchmark:
            "Result/Benchmark/{sample}_Bulkdedup.benchmark" 
        shell:
            "awk  '!a[$1,$2,$3]++' {input.frag_count} > {output.frag_dedup};"


if config["deduplication"] == "cell level":
    rule scatac_fragmentindex:
        input:
            frag = "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv"
        output:
            fraggz = "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv.gz",
            fragindex = "Result/minimap2/{sample}/fragments_corrected_cell_dedup_count.tsv.gz.tbi"
        benchmark:
            "Result/Benchmark/{sample}_FragmentIndex.benchmark"
        shell:
            "bgzip -c {input.frag} > {output.fraggz};"
            "tabix -p bed {output.fraggz}"
elif config["deduplication"] == "bulk level":
    rule scatac_fragmentindex:
        input:
            frag = "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv"
        output:
            fraggz = "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv.gz",
            fragindex = "Result/minimap2/{sample}/fragments_corrected_bulk_dedup_count.tsv.gz.tbi"
        benchmark:
            "Result/Benchmark/{sample}_FragmentIndex.benchmark"
        shell:
            "bgzip -c {input.frag} > {output.fraggz};"
            "tabix -p bed {output.fraggz}" 


