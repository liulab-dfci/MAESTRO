
rule scrna_samplebam:
	input:
		bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
	output:
		bamsample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam",
		baisample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam.bai",
	benchmark:
		"Result/Benchmark/{sample}/{sample}_BamSample.benchmark"
	threads:
		config["cores"]
	shell:
		"""
		samtools view -@ {threads} -s 0.01 -b -o {output.bamsample} {input.bam};
		samtools index -@ {threads} {output.bamsample}
		"""

rule scrna_rseqc_readqual:
	input:
		bamsample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam",
	output:
		qual = "Result/QC/{sample}/{sample}.qual.r"
	params:
		outdirpre = "Result/QC/{sample}/{sample}"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcReadqual.benchmark"
	shell:
		"""
		read_quality.py -i {input.bamsample} -o {params.outdirpre}
		"""

rule scrna_rseqc_nvc:
	input:
		bamsample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam",
	output:
		nvc = "Result/QC/{sample}/{sample}.NVC.xls"
	params:
		outdirpre = "Result/QC/{sample}/{sample}"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcNVC.benchmark"
	shell:
		"""
		read_NVC.py -i {input.bamsample} -o {params.outdirpre}
		"""

rule scrna_rseqc_gc:
	input:
		bamsample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam",
		baisample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam.bai",
	output:
		gc = "Result/QC/{sample}/{sample}.GC.xls"
	params:
		outdirpre = "Result/QC/{sample}/{sample}"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcGC.benchmark"
	shell:
		"""
		read_GC.py -i {input.bamsample} -o {params.outdirpre}
		"""

rule scrna_rseqc_bamstat:
	input:
		bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		bai = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
	output:
		stat = "Result/QC/{sample}/{sample}_bam_stat.txt"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcBamstat.benchmark"
	shell:
		"""
		bam_stat.py -i {input.bam} > {output.stat}
		"""

rule scrna_rseqc_distr:
	input:
		bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
		bai = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
		genome = "%s/annotations/%s_RefSeq.bed" %(SCRIPT_PATH, config["species"]),
	output:
		distr = "Result/QC/{sample}/{sample}_read_distribution.txt"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcDistr.benchmark"
	shell:
		"""
		read_distribution.py -i {input.bam} -r {input.genome} > {output.distr}
		"""

rule scrna_rseqc_genecov:
	input:
		bamsample = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.sample.bam",
		hkgene = "%s/annotations/%s_HouseKeepingGenes.bed" %(SCRIPT_PATH, config["species"]),
	output:
		genecov = "Result/QC/{sample}/{sample}.geneBodyCoverage.txt",
	params:
		outdirpre = "Result/QC/{sample}/{sample}"
	benchmark:
		"Result/Benchmark/{sample}/{sample}_RseqcGenecov.benchmark"
	shell:
		"""
		geneBody_coverage.py -r {input.hkgene} -i {input.bamsample} -o {params.outdirpre}
		"""

rule scrna_rseqc_plot:
	input:
		stat = "Result/QC/{sample}/{sample}_bam_stat.txt",
		distr = "Result/QC/{sample}/{sample}_read_distribution.txt",
		qual = "Result/QC/{sample}/{sample}.qual.r",
		nvc = "Result/QC/{sample}/{sample}.NVC.xls",
		gc = "Result/QC/{sample}/{sample}.GC.xls",
		genecov = "Result/QC/{sample}/{sample}.geneBodyCoverage.txt",
		# countgene = "Result/QC/{sample}/{sample}_count_gene_stat.txt",
	output:
		readdistrplot = "Result/QC/{sample}/{sample}_scRNA_read_distr.png",
		qualplot = "Result/QC/{sample}/{sample}_scRNA_read_quality.png",
		nvcplot = "Result/QC/{sample}/{sample}_scRNA_NVC.png",
		gcplot = "Result/QC/{sample}/{sample}_scRNA_GCcontent.png",
		genecovplot = "Result/QC/{sample}/{sample}_scRNA_genebody_cov.png",
	params:
		outpre = "{sample}/{sample}",
		outdir = "Result/QC",
		rseqc = "TRUE",
		stat = "{sample}/{sample}_bam_stat.txt",
		distr = "{sample}/{sample}_read_distribution.txt",
		qual = "{sample}/{sample}.qual.r",
		nvc = "{sample}/{sample}.NVC.xls",
		gc = "{sample}/{sample}.GC.xls",
		genecov = "{sample}/{sample}.geneBodyCoverage.txt",
		# countgene = "{sample}/{sample}_count_gene_stat.txt",
		# count = config["cutoff"]["count"],
		# gene = config["cutoff"]["gene"]
	benchmark:
		"Result/Benchmark/{sample}/{sample}_QCPlot.benchmark"
	shell:
		"""
		Rscript {RSCRIPT_PATH}/scRNAseq_qc.R --prefix {params.outpre} --outdir {params.outdir} \
		--bamstat {params.stat} --readdistr {params.distr} --qual {params.qual} --nvc {params.nvc} \
		--gc {params.gc} --genecov {params.genecov}
		"""