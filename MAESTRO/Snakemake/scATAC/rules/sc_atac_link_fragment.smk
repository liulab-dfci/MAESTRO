if config["gzip"]:
	rule scatac_link_fragments:
		input:
			frag = lambda wildcards: FILES[wildcards.sample]
		output:
			frag_sort = temp("Result/mapping/{sample}/fragments_sorted_corrected_dedup_count.tsv")
			frag_dedup = "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv",
			fraggz = "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv.gz"
		shell:
			"gunzip -c {input.frag} | sort -k1,1 -k2,2n > {output.frag_sort};"
			"python " + SCRIPT_PATH + "/scATAC_FragmentReshape.py -F {output.frag_sort} -O {output.frag_dedup};"
			"cat {output.frag_dedup} | bgzip > {output.fraggz};"
			"tabix -p bed {output.fraggz}"

else:
	rule scatac_link_fragments:
		input:
			frag = lambda wildcards: FILES[wildcards.sample]
		output:
			frag_sort = temp("Result/mapping/{sample}/fragments_sorted_corrected_dedup_count.tsv")
			frag_dedup = "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv",
			fraggz = "Result/mapping/{sample}/fragments_corrected_dedup_count.tsv.gz"
		shell:
			"sort -k1,1 -k2,2n {input.frag} > {output.frag_sort};"
			"python " + SCRIPT_PATH + "/scATAC_FragmentReshape.py -F {input} -O {output.frag};"
			"cat {output.frag_dedup} | bgzip > {output.fraggz};"
			"tabix -p bed {output.fraggz}"
