rule scatac_map:
    input:
        fasta = config["genome"]["fasta"],
        r1 = lambda wildcards: FILES[wildcards.fastqid]["1"],
        r2 = lambda wildcards: FILES[wildcards.fastqid]["2"],
    output:
        bam = temp("Result/minimap2/{fastqid}.sortedByPos.bam")
    threads:
        config["cores"]
    shell:
        "minimap2 -ax sr -t {threads} {input.fasta} {input.r1} {input.r2} "
        "| samtools view --threads {threads} -b "
        "| samtools sort --threads {threads} -o {output.bam}"

rule scatac_bamrmdp:
    input:
        bam = "Result/minimap2/{fastqid}.sortedByPos.bam"
    output:
        bam = "Result/minimap2/{fastqid}.bam",
        metric = "Result/minimap2/{fastqid}.sortedByPos.rmdp.txt",
        tmp = temp(directory("Result/Tmp/{fastqid}"))
    shell:
        "picard MarkDuplicates INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.metric} TMP_DIR={output.tmp};"
        "rm {input.bam}"

rule scatac_bammerge:
    input:
        bam = expand("Result/minimap2/{fastqid}.bam")
    output:
        bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.bam" ,
        bamlist = "Result/minimap2/sample_bamlist.txt" ,
        fragbed = "Result/QC/sample_frag.bed"
    params:
        sam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.sample.sam" ,
        bamprefix = "Result/minimap2/sample_bamlist_" ,
        subprefix = "Result/minimap2/sample"
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/sample_BamMerge.benchmark"
    shell:
        "ls Result/minimap2/*.bam > {output.bamlist};"
        "split -1000 -d {output.bamlist} {params.bamprefix};"
        "for file in $(ls {params.bamprefix}*); do sub=${{file#{params.bamprefix}}};"
        "samtools merge --threads {threads} -r {params.subprefix}.${{sub}}.sortedByPos.rmdp.bam -b ${{file}};"
        "sinto tagtotag -b {params.subprefix}.${{sub}}.sortedByPos.rmdp.bam --from RG --to CB --delete -o {params.subprefix}.${{sub}}.sortedByPos.rmdp.CBadded.bam -O b; done;"
        "samtools merge --threads {threads} {output.bam} {params.subprefix}.*.sortedByPos.rmdp.CBadded.bam;"
        "rm {params.subprefix}.[0-9]*.sortedByPos.rmdp.bam;"
        "rm {params.subprefix}.[0-9]*.sortedByPos.rmdp.CBadded.bam;"
        "samtools view -@ {threads} -s 0.01 -o {params.sam} {output.bam};"
        "awk '{{if ($9>0) print $9}}' {params.sam} > {output.fragbed};"

rule scatac_fragmentgenerate:
    input:
        bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.bam" ,
    output:
        fragments = "Result/minimap2/fragments_corrected_count.tsv",
    params:
        outdir = "Result/minimap2"
    benchmark:
        "Result/Benchmark/sample_FragGenerate.benchmark"
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentGenerate.py -B {input.bam} -O {params.outdir} --CBtag CB --count"

rule scatac_qcstat_singlecell:
    input:
        bam = "Result/minimap2/{fastqid}.bam",
        promoter = "%s/annotations/%s_promoter.bed" %(SCRIPT_PATH, config["species"]),
        chrM = "%s/annotations/%s_chrM.bed" %(SCRIPT_PATH, config["species"]),
        peak = "Result/Analysis/sample_all_peaks.narrowPeak"
    output:
        log = "Result/Log/bamLog/{fastqid}.mapping.log",
        bam = "Result/minimap2/{fastqid}.sortedByPos.rmdp.unique.bam",
        bed = "Result/minimap2/{fastqid}.sortedByPos.rmdp.unique.bed",
    shell:
        # "samtools flagstat {input.bam} > {output.log};"
        "samtools view -F 2316 -f 0x2 -q 30 -b -o {output.bam} {input.bam};"
        "samtools view {output.bam} -c >> {output.log};"
        "bedtools bamtobed -i {output.bam} > {output.bed};"
        # "grep 'chrM' {output.bed} -c >> {output.log} || true;"
        "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.log} || true;"
        # "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.log} || true ;"

rule scatac_qcstat_bulk:
    input:
        bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.bam" ,
        promoter = SCRIPT_PATH + "/annotations/%s_promoter.bed" %(config["species"]),
        peak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
    output:
        bulk_stat = "Result/QC/flagstat.txt",
        bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.unique.bam" ,
        bed = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.unique.bed" ,
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/sample_BulkQCStat.benchmark"
    shell:
        "samtools flagstat --threads {threads} {input.bam} > {output.bulk_stat};"
        "samtools view -F 2316 -f 0x2 -q 30 -b -o {output.bam} {input.bam};"
        "samtools view {output.bam} -c >> {output.bulk_stat};"
        "bedtools bamtobed -i {output.bam} > {output.bed};"
        "grep 'chrM' {output.bed} -c >> {output.bulk_stat} || true;"
        "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.promoter} -u | wc -l >> {output.bulk_stat} || true;"
        "grep -v 'chrM' {output.bed} | bedtools intersect -wa -a - -b {input.peak} -u | wc -l >> {output.bulk_stat} || true ;"

rule scatac_qcstat_singlemerge:
    input:
        log = expand("Result/Log/bamLog/{fastqid}.mapping.log"),
        # unique = "Result/QC/" + config["outprefix"] + "_uniquereads.txt"
    output:
        stat = "Result/QC/singlecell.txt",
    params:
        log = "Result/Log/bamLog/",
        outdir = "Result/QC"
        # unique = config["outprefix"] + "_uniquereads.txt",
    benchmark:
        "Result/Benchmark/sample_QCMerge.benchmark"
    shell:
        "python " + SCRIPT_PATH + "/scATAC_microfluidic_QC.py --log-dir {params.log} --directory {params.outdir};"

rule scatac_allpeakcall:
    input:
        bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.bam"
    output:
        peak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
        bdg = "Result/Analysis/sample_all_treat_pileup.bdg" ,
    params:
        name = "sample_all" ,
        genome = macs2_genome
    log:
        "Result/Log/sample_macs2_allpeak.log"
    benchmark:
        "Result/Benchmark/sample_AllPeakCall.benchmark"
    shell:
        "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis/ -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.bam}"

if config["shortpeaks"]:
    rule scatac_shortfragment:
        input:
            bam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.bam"
        output:
            shortbam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.150bp.bam"
        threads:
            config["cores"]
        benchmark:
            "Result/Benchmark/sample_ShortFrag.benchmark"
        shell:
            "samtools view -@ {threads} -h {input.bam} | "
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($9)<=150) print}}' | "
            "samtools view -@ {threads} -b -o {output.shortbam}"

    rule scatac_shortpeakcall:
        input:
            shortbam = "Result/minimap2/sample.sortedByPos.rmdp.CBadded.150bp.bam"
        output:
            bed = "Result/Analysis/sample_150bp_peaks.narrowPeak"
        params:
            name = "sample_150bp" ,
            genome = macs2_genome
        log:
            "Result/Log/" + config["outprefix"] + "_macs2_shortpeak.log"
        benchmark:
            "Result/Benchmark/sample_ShortPeakCall.benchmark"
        shell:
            "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.shortbam}"

if config["custompeaks"] and config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
            shortpeak = "Result/Analysis/sample_150bp_peaks.narrowPeak" ,
            custompeak = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/sample_final_peaks.bed"
        params:
            catpeaksort = "Result/Analysis/sample_cat_peaks.bed"
        benchmark:
            "Result/Benchmark/sample_PeakMerge.benchmark"
        shell:
            "cat {input.allpeak} {input.shortpeak} {input.custompeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
elif config["custompeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
            custompeaks = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/sample_final_peaks.bed"
        params:
            catpeaksort = "Result/Analysis/sample_cat_peaks.bed"
        benchmark:
            "Result/Benchmark/sample_PeakMerge.benchmark"
        shell:
            "cat {input.allpeak} {input.custompeaks} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
elif config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
            shortpeak = "Result/Analysis/sample_150bp_peaks.narrowPeak"
        output:
            finalpeak = "Result/Analysis/sample_final_peaks.bed"
        params:
            catpeaksort = "Result/Analysis/sample_cat_peaks.bed"
        benchmark:
            "Result/Benchmark/sample_PeakMerge.benchmark"
        shell:
            "cat {input.allpeak} {input.shortpeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
else:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/sample_all_peaks.narrowPeak" ,
        output:
            finalpeak = "Result/Analysis/sample_final_peaks.bed"
        params:
            catpeaksort = "Result/Analysis/sample_cat_peaks.bed"
        benchmark:
            "Result/Benchmark/sample_PeakMerge.benchmark"
        shell:
            "cat {input.allpeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"

rule scatac_countpeak:
    input:
        finalpeak = "Result/Analysis/sample_final_peaks.bed" ,
        validbarcode = "Result/QC/sample_scATAC_validcells.txt"
    output:
        counts = "Result/Analysis/sample_peak_count.h5"
    threads:
        config["cores"]
    params:
        bamdir = "Result/minimap2",
        species = config["species"],
        outdir = "Result/Analysis",
        outpre = config["outprefix"]
    benchmark:
        "Result/Benchmark/sample_PeakCount.benchmark"
    shell:
        "python " + SCRIPT_PATH + "/scATAC_microfluidic_PeakCount.py --peak {input.finalpeak} --barcode {input.validbarcode} "
        "--bam-dir {params.bamdir} --directory {params.outdir} --outprefix {params.outpre} --cores {threads} --species {params.species}"

rule scatac_qcplot:
    input:
        fragbed = "Result/QC/sample_frag.bed" ,
        single_stat = "Result/QC/singlecell.txt",
    output:
        qcfrag = "Result/QC/" + config["outprefix"] + "_scATAC_fragment_size.png",
        qcfrip = "Result/QC/" + config["outprefix"] + "_scATAC_cell_filtering.png",
        validbarcode = "Result/QC/" + config["outprefix"] + "_scATAC_validcells.txt",
    params:
        outdir = "Result/QC",
        outpre = config["outprefix"],
        fragbed = "sample_frag.bed" ,
        single_stat = "singlecell.txt",
        counts = config["cutoff"]["count"],
        frip = config["cutoff"]["frip"]
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/sample_QCPlot.benchmark"
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_qc.R --fragment {params.fragbed} --singlestat {params.single_stat} "
        "--countcutoff {params.counts} --fripcutoff {params.frip} --prefix {params.outpre} --outdir {params.outdir}"

rule scatac_qcfilter:
    input:
        counts = "Result/Analysis/sample_peak_count.h5" ,
    output:
        filtercount = "Result/QC/sample_filtered_peak_count.h5" ,
    params:
        outdir = "Result/QC",
        outpre = config["outprefix"],
        peak = config["cutoff"]["peak"],
        cell = config["cutoff"]["cell"],
    benchmark:
        "Result/Benchmark/sample_QCFilter.benchmark"
    shell:
        "MAESTRO scatac-qc --format h5 --peakcount {input.counts} --peak-cutoff {params.peak} --cell-cutoff {params.cell} "
        "--directory {params.outdir} --outprefix {params.outpre}"

rule scatac_genescore:
    input:
        filtercount = "Result/QC/sample_filtered_peak_count.h5" ,
        genebed = "sample/annotations/sample_ensembl.bed" %(SCRIPT_PATH, config["species"])
    output:
        genescore = "Result/Analysis/sample_gene_score.h5"
    params:
        genedistance = config["genedistance"],
        species = config["species"],
        outdir = "Result/Analysis",
        outpre = config["outprefix"],
        rpmodel = config["rpmodel"]
    benchmark:
        "Result/Benchmark/sample_GeneScore.benchmark"
    shell:
        "MAESTRO scatac-genescore --format h5 --peakcount {input.filtercount} --species {params.species} --directory {params.outdir} --outprefix {params.outpre} --model {params.rpmodel}"

rule scatac_fragmentindex:
    input:
        frag = "Result/minimap2/fragments_corrected_count.tsv"
    output:
        fraggz = "Result/minimap2/fragments_corrected_count.tsv.gz",
        fragindex = "Result/minimap2/fragments_corrected_count.tsv.gz.tbi"
    benchmark:
        "Result/Benchmark/sample_FragmentIndex.benchmark"
    shell:
        "bgzip -c {input.frag} > {output.fraggz};"
        "tabix -p bed {output.fraggz}"

rule scatac_analysis:
    input:
        filtercount = "Result/QC/sample_filtered_peak_count.h5" ,
        genescore = "Result/Analysis/sample_gene_score.h5" ,
        fraggz = "Result/minimap2/fragments_corrected_count.tsv.gz"
    output:
        specificpeak = "Result/Analysis/sample_DiffPeaks.tsv" ,
        clusterplot = "Result/Analysis/sample_cluster.png" ,
        # annotateplot = "Result/Analysis/sample_annotated.png" ,
        tflist = "Result/Analysis/sample.PredictedTFTop10.txt" ,
        cellcluster = "Result/Analysis/sample_cell_cluster.txt" ,
        # acannotateplot = "Result/Analysis/sample_CistromeTop_annotated.png" ,
        # ms4a1trackplot = "Result/Analysis/sample_MS4A1_genetrack.png" ,
        # cd3dtrackplot = "Result/Analysis/sample_CD3D_genetrack.png" ,
    params:
        outdir = "Result/Analysis",
        genescore = "sample_gene_score.h5" ,
        outpre = config["outprefix"],
        counts = "../QC/sample_filtered_peak_count.h5" ,
        fraggz = "../minimap2/fragments_corrected_count.tsv.gz",
        giggleannotation = config["giggleannotation"],
        species = config["species"],
        signature = config["signature"],
        method = config["method"],
        annotation = config["annotation"],
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/sample_Analysis.benchmark"
    shell:
        "Rscript " + RSCRIPT_PATH + "/scATACseq_pipe.R --peakcount {params.counts} --rpmatrix {params.genescore} "
        "--species {params.species} --prefix {params.outpre} --annotation {params.annotation} --method {params.method} --signature {params.signature} "
        "--gigglelib {params.giggleannotation} --fragment {params.fraggz} --outdir {params.outdir} --thread {threads}"

checkpoint scatac_fragcluster:
    input:
        frag = "Result/minimap2/fragments_corrected_count.tsv",
        cellcluster = "Result/Analysis/sample_cell_cluster.txt" ,
    output:
        directory("Result/Analysis/Cluster/")
    threads:
        config["cores"]
    benchmark:
        "Result/Benchmark/sample_FragCluster.benchmark"
    shell:
        "python " + SCRIPT_PATH + "/scATAC_FragmentSplit.py --frag {input.frag} --cluster {input.cellcluster} --outdir {output}"

rule scatac_clusterpeakcall:
    input:
        frag = "Result/Analysis/Cluster/{cluster}.bed"
    output:
        peak = "Result/Analysis/Cluster/{cluster}_peaks.narrowPeak",
        bdg = "Result/Analysis/Cluster/{cluster}_treat_pileup.bdg",
    params:
        name = "{cluster}",
        outdir = "Result/Analysis/Cluster",
        genome = macs2_genome
    log:
        "Result/Log/{cluster}_macs2_allpeak.log"
    benchmark:
        "Result/Benchmark/{cluster}_AllPeakCall.benchmark"
    shell:
        "macs2 callpeak -f BEDPE -g {params.genome} --outdir {params.outdir} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag}"

rule scatac_report:
    input:
        # bulkqc = "Result/QC/" + config["outprefix"] + "_bam_stat.txt",
        qcfrag = "Result/QC/sample_scATAC_fragment_size.png" ,
        qcfrip = "Result/QC/sample_scATAC_cell_filtering.png" ,
        counts = "Result/QC/sample_filtered_peak_count.h5" ,
        clusterplot = "Result/Analysis/sample_cluster.png" ,
        # annotateplot = "Result/Analysis/sample_annotated.png" ,
        genescore = "Result/Analysis/sample_gene_score.h5" ,
        tflist = "Result/Analysis/sample.PredictedTFTop10.txt" ,
        # acannotateplot = "Result/Analysis/sample_CistromeTop_annotated.png" ,
        # ms4a1trackplot = "Result/Analysis/sample_MS4A1_genetrack.png" ,
        # cd3dtrackplot = "Result/Analysis/sample_CD3D_genetrack.png" ,
    output:
        summaryreport = "Result/sample_scATAC_report.html" ,
    params:
        outpre = "sample",
        inputpath = config["frag"],
        species = config["species"],
        platform = config["platform"],
        inputformat = config["format"],
        outdir = "Result",
    benchmark:
        "Result/Benchmark/sample_Report.benchmark"
    shell:
        # "cp {input.readdistr} {input.qcmap} {input.qcfrag} {input.qcfrip} {input.clusterplot} {input.annotateplot} {output.outdir};"
        "python " + SCRIPT_PATH + "/scATAC_HTMLReport.py --directory {params.outdir} --outprefix {params.outpre} "
        "--input-path {params.inputpath} --species {params.species} --platform {params.platform} --input-format {params.inputformat}"
