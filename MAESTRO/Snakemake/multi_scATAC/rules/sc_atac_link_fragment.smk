
rule scatac_link_fragments:
	input:
		frag = lambda wildcards: FILES[wildcards.sample]
	output: 
		frag_dedup = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv",
		fraggz = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv.gz"
	shell:
		"""
		# make sure it is sorted 
		gunzip -c {input.frag} | sort -k1,1 -k2,2n > {output.frag_dedup}

		cat {output.frag_dedup} | bgzip > {output.fraggz}
		tabix -p bed {output.fraggz}
		"""
