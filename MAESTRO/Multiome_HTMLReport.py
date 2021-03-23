# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-27 17:46:01
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2021-03-23 20:38:12


import os, sys
import snakemake.report

import argparse as ap

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Generate multiome result report. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--rna-fastq-dir", dest = "rna_fastq_dir", type = str, default = "",  
        help = "Directory where RNA fastq files are stored")
    group_input.add_argument("--atac-fastq-dir", dest = "atac_fastq_dir", type = str, default = "",  
        help = "Directory where ATAC fastq files are stored")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_regulator = parser.add_argument_group("Regulator identification arguments")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "MAESTRO", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")
    # group_output.add_argument("--rseqc", dest = "rseqc", action = "store_true", 
    #     help = "Whether or not to run RSeQC. "
    #     "If set, the pipeline will include the RSeQC part and then takes a longer time. "
    #     "By default (not set), the pipeline will skip the RSeQC part.")

    return parser.parse_args()
    

def main():

    SCRIPT_PATH = os.path.dirname(__file__)

    myparser = CommandLineParser()
    directory = myparser.directory
    outpre = myparser.outprefix
    rna_fastqdir = myparser.rna_fastq_dir
    atac_fastqdir = myparser.atac_fastq_dir
    species = myparser.species
    # rseqc = myparser.rseqc

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    # if rseqc:
    #     report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "multiome_template.html")
    #     report_html_temp = open(report_html_tempfile, "r").read()

    #     rna_cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    #     rna_readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_distr.png"%outpre)
    #     rna_readqualplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_quality.png"%outpre)
    #     rna_nvcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_NVC.png"%outpre)
    #     rna_gcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_GCcontent.png"%outpre)
    #     rna_genecovplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_genebody_cov.png"%outpre)
    #     rna_countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)
    #     rna_genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)
    #     rna_geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)
    #     rna_tdtitle = "Cluster-specific regulator identified by %s" %(method)
    #     rna_tdcolname = "log10(%s score)" %(method)

    #     td_list = []
    #     for line in open(cluster_regulator_file,"r").readlines():
    #         if not line.startswith("Cluster"):
    #             items = line.strip().split("\t")
    #             items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
    #             items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
    #             td_list.append(items_str)
    #     td_str = "\n".join(td_list)

    #     #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
    #     report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "readdistr":readdistrplot_link,"readqual":readqualplot_link, "nvc":nvcplot_link, "gc":gcplot_link, "genecov":genecovplot_link, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtabletitle":tdtitle, "regtablecolname":tdcolname, "regtable":td_str}

    # else:
    report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "multiome_template.html")
    report_html_temp = open(report_html_tempfile, "r").read()

    rna_cluster_regulator_file = "Result/RNA/Analysis/%s.PredictedTFTop10.txt"%outpre

    rna_countgeneplot_link = snakemake.report.data_uri_from_file("Result/RNA/QC/%s_scRNA_cell_filtering.png"%outpre)
    rna_genecluster_link = snakemake.report.data_uri_from_file("Result/RNA/Analysis/%s_cluster.png"%outpre)
    rna_geneannotate_link = snakemake.report.data_uri_from_file("Result/RNA/Analysis/%s_annotated.png"%outpre)
    rna_tdtitle = "Cluster-specific regulator identified by LISA"
    rna_tdcolname = "log10(LISA score)"

    rna_td_list = []
    for line in open(rna_cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            rna_td_list.append(items_str)
    rna_td_str = "\n".join(rna_td_list)

    #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],

    # ATAC
    atac_cluster_regulator_file = "Result/ATAC/Analysis/%s.PredictedTFTop10.txt"%outpre

    atac_fragplot_link = snakemake.report.data_uri_from_file("Result/ATAC/QC/%s_scATAC_fragment_size.png"%outpre)
    # mapplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_mapping_summary.png"%outpre)[0]
    atac_fripplot_link = snakemake.report.data_uri_from_file("Result/ATAC/QC/%s_scATAC_cell_filtering.png"%outpre)
    atac_peakcluster_link = snakemake.report.data_uri_from_file("Result/ATAC/Analysis/%s_cluster.png"%outpre)
    if os.path.exists("Result/ATAC/Analysis/%s_annotated.png"%outpre):
        atac_rpannodisplay = "inline"
        atac_rpannotate_link = snakemake.report.data_uri_from_file("Result/ATAC/Analysis/%s_annotated.png"%outpre)
    else:
        atac_rpannodisplay = "none"
        atac_rpannotate_link = ""
    if os.path.exists("Result/ATAC/Analysis/%s_CistromeTop_annotated.png"%outpre):
        atac_caannodisplay = "inline"
        atac_caannotate_link = snakemake.report.data_uri_from_file("Result/ATAC/Analysis/%s_CistromeTop_annotated.png"%outpre)
    else:
        atac_caannodisplay = "none"
        atac_caannotate_link = ""

    if os.path.exists("Result/ATAC/Analysis/%s_MS4A1_genetrack.png"%outpre):
        atac_ms4a1display = "inline"
        atac_ms4a1track_link = snakemake.report.data_uri_from_file("Result/ATAC/Analysis/%s_MS4A1_genetrack.png"%outpre)
    else:
        atac_ms4a1display = "none"
        atac_ms4a1track_link = ""

    if os.path.exists("Result/ATAC/Analysis/%s_CD3D_genetrack.png"%outpre):
        atac_cd3ddisplay = "inline"
        atac_cd3dtrack_link = snakemake.report.data_uri_from_file("Result/ATAC/Analysis/%s_CD3D_genetrack.png"%outpre)
    else:
        atac_cd3ddisplay = "none"
        atac_cd3dtrack_link = ""

    atac_td_list = []
    for line in open(atac_cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            atac_td_list.append(items_str)
    atac_td_str = "\n".join(atac_td_list)

    atac_readdistrplot_link = snakemake.report.data_uri_from_file("Result/ATAC/QC/%s_scATAC_read_distr.png"%outpre)

    joint_filtering_link = snakemake.report.data_uri_from_file("Result/Multiome/%s_multiome_cell_filtering.png"%outpre)
    joint_clustering_link = snakemake.report.data_uri_from_file("Result/Multiome/%s_cluster_wsnn.png"%outpre)
    joint_annotate_link = snakemake.report.data_uri_from_file("Result/Multiome/%s_annotated_wsnn.png"%outpre)


    report_html_instance = report_html_temp % {"outprefix": outpre, "rna_fastqdir": rna_fastqdir, "atac_fastqdir": atac_fastqdir, "species":species,
    "cellfilter_joint": joint_filtering_link,"cluster_joint":joint_clustering_link, "annotate_joint":joint_annotate_link,
    "countgene_rna":rna_countgeneplot_link, "genecluster_rna":rna_genecluster_link, 
    "geneannotate_rna":rna_geneannotate_link, "regtabletitle_rna":rna_tdtitle, 
    "regtablecolname_rna":rna_tdcolname, "regtable_rna":rna_td_str,
    "readdistr_atac":atac_readdistrplot_link, "fragment_atac":atac_fragplot_link, "frip_atac":atac_fripplot_link, 
    "peakcluster_atac":atac_peakcluster_link, "rpannodisplay_atac": atac_rpannodisplay, "rpannotate_atac":atac_rpannotate_link, "regtable_atac":atac_td_str,
    "caannodisplay_atac": atac_caannodisplay, "caannotate_atac": atac_caannotate_link, 'cd3ddisplay_atac': atac_cd3ddisplay, 'ms4a1display_atac': atac_ms4a1display,
    'cd3dtrack_atac': atac_cd3dtrack_link, "ms4a1track_atac": atac_ms4a1track_link}



    report_html_instancefile = os.path.join(directory, outpre + "_multiome_report.html")
    outf = open(report_html_instancefile,"w")
    outf.write(report_html_instance)
    outf.close()


if __name__ == "__main__":
    main()

