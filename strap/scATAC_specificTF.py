import os
import sys
import argparse
import pandas as pd 


def giggleSearchBed(bed_path, output, assembly):
    # if not os.path.exists(output):
    #     os.makedirs(output)
    prefix = "%s/%s" % (output, bed_path.split("/")[-1])

    ### NOTE: This table should be replaced after we change the way of fetching 5foldReakRP
    ant = pd.read_csv("/project/dev/changxin/strap/gindex_201905/DC_haveProcessed_20190506_filepath.xls", header=0, index_col=False, sep="\t", engine="c")
    ant.index = list(map(str, ant.loc[:, "DCid"].tolist()))
    ant = ant.loc[:, ["factor", "cellLine", "CellType", "Tissue", "5foldPeakRP"]]

    cmd = "sort --buffer-size 2G -k1,1 -k2,2n -k3,3n %s | bgzip -c > %s.gz"%(bed_path, prefix)
    os.system(cmd)

    ### NOTE: Change the path of giggle and giggle index
    os.system("/project/dev/changxin/Software/Giggle/giggle/bin/giggle search -i /project/dev/changxin/project_01/changxin/gindex/%stf_gindex_1k -q %s.gz -s > %s.result.xls"%(assembly, prefix, prefix))
    result_df = pd.read_csv("%s.result.xls"%prefix, sep="\t", index_col=False)
    result_df.index = [i.replace("_5foldPeak.bed.gz", "").split("/")[-1] for i in result_df["#file"]]
    result_df = result_df.loc[:, ["file_size", "overlaps", "combo_score"]]

    os.system("rm %s.gz"%prefix)
    os.system("rm %s.result.xls"%prefix)
    res_df = pd.concat([result_df, ant], axis=1, join="inner").sort_values(by="combo_score", ascending=False)
    res_df.columns = ["sample_peak_number", "overlap_peak_number", "giggle_score", "factor", "cell_line", "cell_type", "tissue", "5foldPeakRP"]
    res_df = res_df.loc[res_df["overlap_peak_number"]>0, :]
    res_df["biological_resource"] = list(map(lambda x: "%s;%s;%s"%(res_df["cell_line"][x], res_df["cell_type"][x], res_df["tissue"][x]), range(res_df.shape[0])))
    res_df["sample_id"] = res_df.index.tolist()
    res_df["species"] = assembly
    res_df = res_df.loc[:, ["sample_id", "factor", "species", "biological_resource", "giggle_score", "sample_peak_number", "overlap_peak_number", "5foldPeakRP"]]
    res_df = res_df.drop_duplicates(subset="factor").iloc[:10, :]
    res_df.index = range(res_df.shape[0])
    res_df.to_csv("%s/giggleTF.txt"%output, sep="\t", header=True, index=False)
    for i in res_df.index:
        df_5fold = pd.read_csv(res_df.loc[i, "5foldPeakRP"], comment="#", sep="\t", header=None, index_col=False)
        df_5fold.columns = ["chrom", "txStart", "txEnd", "refseq", "score","strand","symbol"]
        df_5fold = df_5fold.loc[df_5fold["score"]>0, :]
        df_5fold = df_5fold.sort_values(by="score", ascending=False)
        target_genes = df_5fold.drop_duplicates(subset="symbol").iloc[:500, 6].tolist()
        f = open("%s/%s_%s_target_genes_top500.txt" % (output, res_df.loc[i, "factor"], res_df.loc[i, "sample_id"]), "w")
        f.write("\n".join(target_genes))
        f.close()
    return res_df

def searchClusterSpecificTFs(peak_cluster, outpath, assembly, spt="-"):
    df = pd.read_csv(peak_cluster, sep="\t", header=0, index_col=0)
    bed_files = []
    for cluster in set(df["cluster"]):
        clusteri_pos = df.loc[df["cluster"]==cluster, "gene"].tolist()
        clusteri_str = "\n".join(list(map(lambda x: "\t".join(x.strip().split(spt)), clusteri_pos)))
        output = "%s/cluster_%s/"%(outpath, cluster)
        if not os.path.exists(output):
            os.makedirs(output)
        f = open("%s/cluster_specific_peaks.bed" % output, "w")
        f.write(clusteri_str)
        f.close()
        bed_files.append(output)

    for output in bed_files:
        res_giggle_bed = giggleSearchBed("%s/cluster_specific_peaks.bed" % output, output, assembly)
        print(output)

def main():
### NOTE: Add some print out information for the script
    try:
        parser = argparse.ArgumentParser(description="""Search cluster specific TFs according to scATAC-seq data""")
        parser.add_argument( '-a', '--assembly', dest='assembly', type=str, required=True, choices=['hg38', 'mm10'], help='select either hg38 or mm10 for genome assembly')
        parser.add_argument( '-p', '--peaks', dest='peaks', type=str, required=True, help='cluster specific peaks output from Seurate2 or Seurate3 delimitated by tab')
        parser.add_argument( '-o', '--output', dest='outpath', type=str, required=True, help='directory of output file')
        parser.add_argument( '-s', '--separate', dest='separate', required=True, help='_ for Seurat2 and - for Seurat3')
        args = parser.parse_args()
        searchClusterSpecificTFs(args.peaks, args.outpath, args.assembly, args.separate)
    except KeyboardInterrupt:
        sys.stderr.write("User interrupted me!\n")
        sys.exit(0)

if __name__ == '__main__':
    main()
