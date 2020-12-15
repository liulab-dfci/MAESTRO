

if config["clusterpeaks"]:
    checkpoint scatac_fragcluster:
        input: 
            merged_frag = "Result/Analysis/Batch/all_samples_fragments_corrected_dedup_count.tsv.gz",
            cellcluster = "Result/Analysis/Batch/all_samples_cell_cluster.txt"
        output:
            directory("Result/Analysis/Cluster/per_cluster")
        log:
            "Result/Log/scatac_fragcluster.log"
        benchmark:
            "Result/Benchmark/all_sample_split_frag_by_cluster.benchmark"

        shell:
            "python utils/scATAC_split_by_cluster.py -F {input.merged_frag} -C {input.cellcluster} -S by_cluster -O {output} 2> {log}"


    checkpoint scatac_fragcluster_per_sample:
        input: 
            merged_frag = "Result/Analysis/Batch/all_samples_fragments_corrected_dedup_count.tsv.gz",
            cellcluster = "Result/Analysis/Batch/all_samples_cell_cluster.txt"
        output:
            directory("Result/Analysis/Cluster/per_cluster_sample")
        log:
            "Result/Log/scatac_fragcluster_per_sample.log"    
        benchmark:
            "Result/Benchmark/all_sample_split_frag_by_cluster.benchmark"

        shell:
            "python utils/scATAC_split_by_cluster.py -F {input.merged_frag} -C {input.cellcluster} -S by_sample_cluster -O {output} 2> {log}"

    rule scatac_call_peak_per_cluster:
        input:
            "Result/Analysis/Cluster/per_cluster/{cluster}.tsv"
        output:
            bdg = "Result/Analysis/Cluster/per_cluster/{cluster}_treat_pileup.bdg"
        params:
            name = "{cluster}",
            genome = macs2_genome
        shell:
            "macs2 callpeak -f BEDPE -g {params.genome} --outdir Result/Analysis/Cluster/per_cluster/ -n {params.name} -B -q 0.05 --nomodel --extsize=50 --SPMR --keep-dup all -t {input}"

    rule scatac_bdg2bw_per_cluster:
        input: 
            bdg = "Result/Analysis/Cluster/per_cluster/{cluster}_treat_pileup.bdg"
        output:
            bw = "Result/Analysis/Cluster/per_cluster/{cluster}.bw"
        params:
            chrom_len = config["chrom_len"],
            sort_bdg = "Result/Analysis/Cluster/per_cluster/{cluster}_sort.bdg",
            clip_bdg = "Result/Analysis/Cluster/per_cluster/{cluster}_sort.clip"
        shell:
            """
            # https://gist.github.com/taoliu/2469050  bdg file generated from MACS2 can exceed chromosome limits
            bedtools slop -i {input.bdg} -g {params.chrom_len} -b 0 | ./utils/bedClip stdin {params.chrom_len} {params.clip_bdg}
            sort -k1,1 -k2,2n {params.clip_bdg} > {params.sort_bdg}

            ./utils/bedGraphToBigWig {params.sort_bdg} {params.chrom_len} {output.bw}

            rm {params.sort_bdg}
            rm {params.clip_bdg}
            """

    ## make RPM normalized bw per cluster per sample without using MACS 
    rule scatac_frag2bw_per_cluster_sample:
        input: 
            tsv = "Result/Analysis/Cluster/per_cluster_sample/{cluster}@{sample}.tsv"
        output:
            bw = "Result/Analysis/Cluster/per_cluster_sample/{cluster}@{sample}.bw"
        params:
            chrom_len = config["chrom_len"],
            bdg = "Result/Analysis/Cluster/per_cluster_sample/{cluster}@{sample}.bdg"
        shell:
            """
            total=`wc -l {input.tsv} | awk '{{print $1}}'`
            # make RPM normalized bigwig
            bedtools genomecov -i {input.tsv} -bg -g {params.chrom_len} | sort -k1,1 -k2,2n - | awk -v OFS='\t' -v total=$total '{{print $1,$2,$3,$4*1000000/total}}' > {params.bdg}

            ./utils/bedGraphToBigWig {params.bdg} {params.chrom_len} {output.bw}

            rm {params.bdg}
            """





