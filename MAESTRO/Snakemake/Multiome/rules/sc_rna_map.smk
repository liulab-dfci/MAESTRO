# rule for sequencing mapping

from MAESTRO.scRNA_utility import getfastq_10x

rule scrna_map:
    input:
        mapindex = config["genome"]["rna_mapindex"],
        whitelist = config["barcode"]["rna_whitelist"],
        fastqs = config["rna_fastqdir"],
    output:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
        bai = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam.bai" %(config["outprefix"]),
        rawmtx = "Result/RNA/STAR/%sSolo.out/Gene/raw/matrix.mtx" %(config["outprefix"]),
        feature = "Result/RNA/STAR/%sSolo.out/Gene/raw/features.tsv" %(config["outprefix"]),
        barcode = "Result/RNA/STAR/%sSolo.out/Gene/raw/barcodes.tsv" %(config["outprefix"]),
    params:
        outdir = "Result/RNA/STAR/" + config["outprefix"],
        transcript = getfastq_10x(config["rna_fastqdir"], config["rna_fastqprefix"])["transcript"],
        barcode = getfastq_10x(config["rna_fastqdir"], config["rna_fastqprefix"])["barcode"],
        decompress = getfastq_10x(config["rna_fastqdir"], config["rna_fastqprefix"])["decompress"],
        barcodestart = 1,
        barcodelength = 16,
        umistart = 17,
        umilength = 12,
    log:
        "Result/Log/%s_scRNA_STAR.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_scRNA_STAR.benchmark" %(config["outprefix"])
    threads:
        config["cores"]
    shell:
        "STAR --runMode alignReads --genomeDir {input.mapindex} --runThreadN {threads} "
        "--outFileNamePrefix {params.outdir} --outSAMtype BAM SortedByCoordinate "
        "--soloType CB_UMI_Simple --soloCBwhitelist {input.whitelist} "
        "--soloCBstart {params.barcodestart} --soloCBlen {params.barcodelength} "
        "--soloUMIstart {params.umistart} --soloUMIlen {params.umilength} "
        "--soloCBmatchWLtype 1MM_multi_pseudocounts --soloUMIfiltering MultiGeneUMI "
        "--readFilesIn {params.transcript} {params.barcode} --readFilesCommand {params.decompress} "
        "> {log};"
        "samtools index -b -@ {threads} {output.bam}"
