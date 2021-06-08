

STAR_VERSION = subprocess.check_output("STAR --version", shell=True)

if config["platform"] == "10x-genomics":
    rule STARsolo:
        """
        STARsolo to align single cell or single nuceli data
        specify --soloFeatures Gene for single-cell  data
        specify --soloFeatures GeneFull for single-nuclei data
        specify --soloFeatures Gene GeneFull for getting both counts in exons level and exon + intron level (velocity)
        """
        input: 
            mapindex = config["genome"]["mapindex"],
            #gtf = config["genome"]["gtf"],
            whitelist = config["barcode"]["whitelist"]
        output: 
            bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            bai = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
            rawmtx = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/matrix.mtx",
            feature = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/features.tsv",
            barcode = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/barcodes.tsv"
        params:
            star_custom = config.get("STARsolo_custom", ""),
            outprefix = "Result/STAR/{sample}/{sample}",
            transcript = lambda wildcards: ','.join(FILES[wildcards.sample]["R2"]),
            barcode = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
            barcodestart = config["barcode"]["barcodestart"],
            barcodelength = config["barcode"]["barcodelength"],
            umistart = config["barcode"]["umistart"],
            umilength = config["barcode"]["umilength"]
        version: STAR_VERSION        
        log:
            "Result/Log/{sample}_STAR.log" 
        benchmark:
            "Result/Benchmark/{sample}_STAR.benchmark"
        threads:
            config.get("STARsolo_threads", "")

        shell:
            """
            
            STAR \
                --runMode alignReads \
                --genomeDir {input.mapindex} \
                --runThreadN {threads} \
                --outFileNamePrefix {params.outprefix} \
                --outSAMtype BAM SortedByCoordinate \
                --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
                --soloType CB_UMI_Simple \
                {params.star_custom} \
                --soloCBwhitelist {input.whitelist} \
                --soloCBstart {params.barcodestart} \
                --soloCBlen {params.barcodelength} \
                --soloUMIstart {params.umistart} \
                --soloUMIlen {params.umilength} \
                --soloCBmatchWLtype 1MM_multi_pseudocounts \
                --soloUMIfiltering MultiGeneUMI \
                --readFilesIn {params.transcript} {params.barcode} \
				--readFilesCommand zcat \
                --genomeSAindexNbases 2 \
                > {log} 2>&1
                
            samtools index -b -@ {threads} {output.bam} >> {log} 2>&1
            """
elif config["platform"] == "Dropseq": 
    rule scrna_map:
        input:
            mapindex = config["genome"]["mapindex"],
            whitelist = config["barcode"]["whitelist"]
        output:
            bam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            bai = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
            rawmtx = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/matrix.mtx",
            feature = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/features.tsv",
            barcode = "Result/STAR/{sample}/{sample}Solo.out/Gene/raw/barcodes.tsv"
        params:
            outprefix = "Result/STAR/{sample}/{sample}",
            transcript = lambda wildcards: ','.join(FILES[wildcards.sample]["R2"]),
            barcode = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
            decompress = getfastq_dropseq(config["fastqdir"], config["fastq"]["barcode"], config["fastq"]["transcript"])["decompress"],
            barcodestart = config["barcode"]["barcodestart"],
            barcodelength = config["barcode"]["barcodelength"],
            umistart = config["barcode"]["umistart"],
            umilength = config["barcode"]["umilength"]
        version: STAR_VERSION
        log:
            "Result/Log/%s_STAR.log" %(config["outprefix"])
        threads:
            config["cores"]
        benchmark:
            "Result/Benchmark/%s_STAR.benchmark" %(config["outprefix"])
        shell:
            """
            STAR \
                --runMode alignReads \
                --genomeDir {input.mapindex} \
                --runThreadN {threads} \
                --outFileNamePrefix {params.outprefix} \
                --outSAMtype BAM SortedByCoordinate \
                --soloType CB_UMI_Simple \
                --soloCBwhitelist {input.whitelist} \
                --soloCBstart {params.barcodestart} \
                --soloCBlen {params.barcodelength} \
                --soloUMIstart {params.umistart} \
                --soloUMIlen {params.umilength} \
                --soloCBmatchWLtype 1MM_multi_pseudocounts \
                --soloUMIfiltering MultiGeneUMI \
                --readFilesIn {params.transcript} {params.barcode} \
                --readFilesCommand {params.decompress} \
                > {log} 2>&1
            
            samtools index -b -@ {threads} {output.bam}
            """
elif config["platform"] == "Smartseq2":
    rule scrna_map:
        input:
            mapindex = config["genome"]["mapindex"],
            fastq1 = lambda wildcards: ','.join(FILES[wildcards.sample]["R1"]),
            fastq2 = lambda wildcards: ','.join(FILES[wildcards.sample]["R2"])
        output:
            genomebam = "Result/STAR/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            transbam = "Result/STAR/{sample}/{sample}Aligned.toTranscriptome.out.bam"
        params:
            outdir = "Result/STAR/{sample}"
        version: STAR_VERSION
        log:
            "Result/Log/{sample}/{sample}_STAR_map.log"
        threads:
            config["cores"]
        shell:
            """
            STAR \
				--genomeDir {input.mapindex} \
				--runThreadN {threads} \
				--outFilterMultimapNmax 500 \
				--outFilterMismatchNmax 3 \
				--quantMode TranscriptomeSAM \
				--outFileNamePrefix {params.outdir} \
				--outSAMtype BAM SortedByCoordinate \
				--readFilesIn {input.fastq1} {input.fastq2} \
				> {log} 2>&1
            """







