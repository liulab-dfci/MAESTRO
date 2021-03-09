# rule for pre-process scATAC fastq

rule scatac_preprocess:
    input:
        fastqs = config["fastqdir"],
        fasta = config["genome"]["fasta"]
    output:
        r1cat = temp(os.path.join("Result/ATAC/Tmp", "%s_R1.fastq" %(config["fastqprefix"]))),
        r2cat = temp(os.path.join("Result/ATAC/Tmp", "%s_R2.fastq" %(config["fastqprefix"]))),
        r3cat = temp(os.path.join("Result/ATAC/Tmp", "%s_R3.fastq" %(config["fastqprefix"]))),
    params:
        r1 = getfastq_10x(config["fastqdir"], config["fastqprefix"])["r1"],
        r2 = getfastq_10x(config["fastqdir"], config["fastqprefix"])["barcode"],
        r3 = getfastq_10x(config["fastqdir"], config["fastqprefix"])["r3"],
        cmd = getfastq_10x(config["fastqdir"], config["fastqprefix"])["command"],
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_Preprocess.benchmark" %(config["outprefix"])
    shell:
        "{params.cmd} {params.r1} > {output.r1cat};"
        "{params.cmd} {params.r2} > {output.r2cat};"
        "{params.cmd} {params.r3} > {output.r3cat};"

rule scatac_fqaddbarcode:
    input:
        r1 = os.path.join("Result/ATAC/Tmp", "%s_R1.fastq" %(config["fastqprefix"])),
        r2 = os.path.join("Result/ATAC/Tmp", "%s_R2.fastq" %(config["fastqprefix"])),
        r3 = os.path.join("Result/ATAC/Tmp", "%s_R3.fastq" %(config["fastqprefix"])),
    output:
        r1 = temp(os.path.join("Result/ATAC/Tmp", "%s_R1.barcoded.fastq" %(config["fastqprefix"]))),
        r3 = temp(os.path.join("Result/ATAC/Tmp", "%s_R3.barcoded.fastq" %(config["fastqprefix"]))),
        # r1 = "%s/%s_R1.barcoded.fastq" %(config["fastqdir"], config["fastqprefix"]),
        # r3 = "%s/%s_R3.barcoded.fastq" %(config["fastqdir"], config["fastqprefix"]),
    benchmark:
        "Result/Benchmark/%s_scATAC_FqAddbarcode.benchmark" %(config["outprefix"])
    shell:
        "base=`head -n 2 {input.r2} | tail -n 1 | wc -L`;"
        "sinto barcode --barcode_fastq {input.r2} --read1 {input.r1} --read2 {input.r3} -b $base;"

