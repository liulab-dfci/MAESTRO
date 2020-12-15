_shortfragment_threads = 2

macs2_genome = "hs" if config["species"] == "GRCh38" else "mm"

rule scatac_allpeakcall:
    input:
        frag_dedup = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv" 
    output:
        peak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
        bdg = "Result/Analysis/{sample}/{sample}_all_treat_pileup.bdg"
    params:
        name = "{sample}_all",
        genome = macs2_genome
    log:
        "Result/Log/{sample}_macs2_allpeak.log" 
    benchmark:
        "Result/Benchmark/{sample}_AllPeakCall.benchmark" 
    shell:
        "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag_dedup}"

if config["shortpeaks"]:
    rule scatac_shortfragment:
        input:
            frag_dedup = "Result/minimap2/{sample}/fragments_corrected_dedup_count.tsv"
        output:
            frag_short = "Result/minimap2/{sample}/fragments_corrected_150bp.tsv" 
        threads:
            _shortfragment_threads
        benchmark:
            "Result/Benchmark/{sample}_ShortFrag.benchmark" 
        shell:
            "awk -F'\\t' 'function abs(x){{return ((x < 0.0) ? -x : x)}} {{if (abs($3-$2)<=150) print}}' {input.frag_dedup} > {output.frag_short}" 
    
    rule scatac_shortpeakcall:
        input:
            frag_short = "Result/minimap2/{sample}/fragments_corrected_150bp.tsv" 
        output:
            bed = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak"
        params:
            name = "{sample}_150bp",
            genome = macs2_genome 
        log:
            "Result/Log/{sample}_macs2_shortpeak.log" 
        benchmark:
            "Result/Benchmark/{sample}_ShortPeakCall.benchmark" 
        shell:
            "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/Analysis/{wildcards.sample} -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input.frag_short}"

if config["custompeaks"] and config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak" ,
            shortpeak = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak",
            custompeak = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} {input.shortpeak} {input.custompeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["custompeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
            custompeaks = config["custompeaksloc"]
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed"  
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark"          
        shell:
            """
            cat {input.allpeak} {input.custompeaks} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

elif config["shortpeaks"]:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak",
            shortpeak = "Result/Analysis/{sample}/{sample}_150bp_peaks.narrowPeak" 
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} {input.shortpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

else:
    rule scatac_mergepeak:
        input:
            allpeak = "Result/Analysis/{sample}/{sample}_all_peaks.narrowPeak" 
        output:
            finalpeak = "Result/Analysis/{sample}/{sample}_final_peaks.bed" 
        params:
            catpeaksort = "Result/Analysis/{sample}/{sample}_cat_peaks.bed" 
        benchmark:
            "Result/Benchmark/{sample}_PeakMerge.benchmark" 
        shell:
            """
            cat {input.allpeak} \
            | sort -k1,1 -k2,2n | cut -f 1-4 > {params.catpeaksort}

            mergeBed -i {params.catpeaksort} | grep -v '_' | grep -v 'chrEBV' > {output.finalpeak}

            rm {params.catpeaksort}
            """

rule scatac_bdg2bw:
    input: 
        bdg = "Result/Analysis/{sample}/{sample}_all_treat_pileup.bdg"
    output:
        bw = "Result/Analysis/{sample}/{sample}.bw"
    benchmark:
        "Result/Benchmark/{sample}_Bdg2Bw.benchmark"
    params:
        chrom_len = config["chrom_len"],
        sort_bdg = "Result/Analysis/{sample}/{sample}_sort.bdg",
        clip_bdg = "Result/Analysis/{sample}/{sample}_sort.clip"
    shell:
        """
        # https://gist.github.com/taoliu/2469050  bdg file generated from MACS2 can exceed chromosome limits
        bedtools slop -i {input.bdg} -g {params.chrom_len} -b 0 | ./utils/bedClip stdin {params.chrom_len} {params.clip_bdg}
        sort -k1,1 -k2,2n {params.clip_bdg} > {params.sort_bdg}

        ./utils/bedGraphToBigWig {params.sort_bdg} {params.chrom_len} {output.bw}

        rm {params.sort_bdg}
        rm {params.clip_bdg}
        """



