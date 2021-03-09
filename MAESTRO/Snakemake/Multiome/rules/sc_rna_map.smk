# rule for sequencing mapping

from MAESTRO.scRNA_utility import getfastq_10x

rule scrna_map:
    input:
        mapindex = config["genome"]["mapindex"],
        whitelist = config["barcode"]["whitelist"],
        fastqs = config["fastqdir"],
    output:
        bam = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam" %(config["outprefix"]),
        bai = "Result/RNA/STAR/%sAligned.sortedByCoord.out.bam.bai" %(config["outprefix"]),
        rawmtx = "Result/RNA/STAR/%sSolo.out/Gene/raw/matrix.mtx" %(config["outprefix"]),
        feature = "Result/RNA/STAR/%sSolo.out/Gene/raw/features.tsv" %(config["outprefix"]),
        barcode = "Result/RNA/STAR/%sSolo.out/Gene/raw/barcodes.tsv" %(config["outprefix"]),
    params:
        outdir = "Result/RNA/STAR/" + config["outprefix"],
        transcript = getfastq_10x(config["fastqdir"], config["fastqprefix"])["transcript"],
        barcode = getfastq_10x(config["fastqdir"], config["fastqprefix"])["barcode"],
        decompress = getfastq_10x(config["fastqdir"], config["fastqprefix"])["decompress"],
        barcodestart = config["barcode"]["barcodestart"],
        barcodelength = config["barcode"]["barcodelength"],
        umistart = config["barcode"]["umistart"],
        umilength = config["barcode"]["umilength"],
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
