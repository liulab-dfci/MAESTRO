# rule for peak calling

macs2_genome = "hs" if config["species"] == "GRCh38" else "mm"

rule scatac_allpeakcall:
    input:
        bam = "Result/ATAC/minimap2/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"])
    output:
        peak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
        bdg = "Result/ATAC/Analysis/%s_all_treat_pileup.bdg" %(config["outprefix"]),
    params:
        name = "%s_all" %(config["outprefix"]),
        genome = macs2_genome,
        outdir = "Result/ATAC/Analysis"
    log:
        "Result/Log/%s_scATAC_macs2_allpeak.log" %(config["outprefix"])
    benchmark:
        "Result/Benchmark/%s_AllPeakCall.benchmark" %(config["outprefix"])
    shell:
        "macs2 callpeak -f BAMPE -g {params.genome} --outdir {params.outdir} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.bam}"


if config["shortpeaks"]:
    rule scatac_shortfragment:
        input:
            bam = "Result/ATAC/minimap2/%s.sortedByPos.rmdp.CBadded.bam" %(config["outprefix"])
        output:
            shortbam = "Result/ATAC/minimap2/%s.sortedByPos.rmdp.CBadded.150bp.bam" %(config["outprefix"])
        threads:
            config["cores"]
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortFrag.benchmark" %(config["outprefix"])
        shell:
            "samtools view -@ {threads} -h {input.bam} | "
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($9)<=150) print}}' | "
            "samtools view -@ {threads} -b -o {output.shortbam}"

    rule scatac_shortpeakcall:
        input:
            shortbam = "Result/ATAC/minimap2/%s.sortedByPos.rmdp.CBadded.150bp.bam" %(config["outprefix"])
        output:
            bed = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"])
        params:
            name = "%s_150bp" %(config["outprefix"]),
            genome = macs2_genome
        log:
            "Result/Log/%s_scATAC_macs2_shortpeak.log" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_ShortPeakCall.benchmark" %(config["outprefix"])
        shell:
            "macs2 callpeak -f BAMPE -g {params.genome} --outdir Result/Analysis -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.shortbam}"


if config["custompeaks"] and config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            shortpeak = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"]),
            custompeak = config["custompeaksloc"]
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            "cat {input.allpeak} {input.shortpeak} {input.custompeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
elif config["custompeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            custompeaks = config["custompeaksloc"]
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            "cat {input.allpeak} {input.custompeaks} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
elif config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
            shortpeak = "Result/ATAC/Analysis/%s_150bp_peaks.narrowPeak" %(config["outprefix"])
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            "cat {input.allpeak} {input.shortpeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
else:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/ATAC/Analysis/%s_all_peaks.narrowPeak" %(config["outprefix"]),
        output:
            finalpeak = "Result/ATAC/Analysis/%s_final_peaks.bed" %(config["outprefix"])
        params:
            catpeaksort = "Result/ATAC/Analysis/%s_cat_peaks.bed" %(config["outprefix"])
        benchmark:
            "Result/Benchmark/%s_scATAC_PeakMerge.benchmark" %(config["outprefix"])
        shell:
            "cat {input.allpeak} "
            "| sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort};"
            "mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak};"
            "rm {params.catpeaksort}"
