# rule for scATAC part


rule scatac_genescore:
    input:
        filtercount = "Result/Multiome/%s_multiome_peak_count.h5" %(config["outprefix"]),
        genebed = "%s/annotations/%s_ensembl.bed" %(SCRIPT_PATH, config["species"])
    output:
        genescore = "Result/ATAC/Analysis/%s_gene_score.h5" %(config["outprefix"])
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/ATAC/Analysis",
        outpre = config["outprefix"],
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/%s_scATAC_GeneScore.benchmark" %(config["outprefix"])
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"

rule scatac_fragmentindex:
    input:
        frag = "Result/ATAC/minimap2/fragments_corrected_count.tsv"
    output:
        fraggz = "Result/ATAC/minimap2/fragments_corrected_count.tsv.gz",
        fragindex = "Result/ATAC/minimap2/fragments_corrected_count.tsv.gz.tbi"
    benchmark:
        "Result/Benchmark/%s_scATAC_FragmentIndex.benchmark" %(config["outprefix"])
    shell:
        "bgzip -c {input.frag} > {output.fraggz};"
        "tabix -p bed {output.fraggz}"

rule scatac_analysis:
    input:
        filtercount = "Result/Multiome/%s_multiome_peak_count.h5" %(config["outprefix"]),
        genescore = "Result/ATAC/Analysis/%s_gene_score.h5" %(config["outprefix"]),
        fraggz = "Result/ATAC/minimap2/fragments_corrected_count.tsv.gz"
    output:
        specificpeak = "Result/ATAC/Analysis/%s_DiffPeaks.tsv" %(config["outprefix"]),
        clusterplot = "Result/ATAC/Analysis/%s_cluster.png" %(config["outprefix"]),
        tflist = "Result/ATAC/Analysis/%s.PredictedTFTop10.txt" %(config["outprefix"]),
        cellcluster = "Result/ATAC/Analysis/%s_cell_cluster.txt" %(config["outprefix"]),
        atacobject = "Result/ATAC/Analysis/%s_scATAC_Object.rds" %(config["outprefix"]),
    params:
        outdir = "Result/ATAC/Analysis",
        genescore = "%s_gene_score.h5" %(config["outprefix"]),
        outpre = config["outprefix"],
        counts = "../../Multiome/%s_multiome_peak_count.h5" %(config["outprefix"]),
        fraggz = "../minimap2/fragments_corrected_count.tsv.gz",
        giggleannotation = config["giggleannotation"],
        species = config["species"],
        signature = config["signature"],
        method = config["method"],
        annotation = config["annotation"],
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_Analysis.benchmark" %(config["outprefix"])
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_pipe.R --peakcount {params.counts} --rpmatrix {params.genescore} "
        "--species {params.species} --prefix {params.outpre} --annotation {params.annotation} --method {params.method} --signature {params.signature} "
        "--gigglelib {params.giggleannotation} --fragment {params.fraggz} --outdir {params.outdir} --thread {threads}"

checkpoint scatac_fragcluster:
    input:
        frag = "Result/ATAC/minimap2/fragments_corrected_count.tsv",
        cellcluster = "Result/ATAC/Analysis/%s_cell_cluster.txt" %(config["outprefix"]),
    output:
        directory("Result/ATAC/Analysis/Cluster/")
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/%s_scATAC_FragCluster.benchmark" %(config["outprefix"])
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentSplit.py --frag {input.frag} --cluster {input.cellcluster} --outdir {output}"

rule scatac_clusterpeakcall:
    input:
        frag = "Result/ATAC/Analysis/Cluster/{cluster}.bed"
    output:
        peak = "Result/ATAC/Analysis/Cluster/{cluster}_peaks.narrowPeak",
        bdg = "Result/ATAC/Analysis/Cluster/{cluster}_treat_pileup.bdg",
    params:
        name = "{cluster}",
        outdir = "Result/ATAC/Analysis/Cluster",
        genome = macs2_genome
    log:
        "Result/Log/{cluster}_scATAC_macs2_allpeak.log"
    benchmark:
        "Result/Benchmark/{cluster}_scATAC_AllPeakCall.benchmark"
    shell:
        "macs2 callpeak -f BEDPE -g {params.genome} --outdir {params.outdir} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag}"
