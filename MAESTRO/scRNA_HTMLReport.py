# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-27 17:46:01
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-29 04:11:13


import os, sys
import snakemake.report

import argparse as ap

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Generate scRNA result report. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--platform", dest = "platform", default = "10x-genomics", 
        choices = ["10x-genomics", "Dropseq", "Smartseq2"], 
        help = "Platform of single cell RNA-seq. DEFAULT: 10x-genomics.")
    group_input.add_argument("--fastq-dir", dest = "fastq_dir", type = str, default = "",  
        help = "Directory where fastq files are stored")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_regulator = parser.add_argument_group("Regulator identification arguments")
    group_regulator.add_argument("--method", dest = "method", type = str, 
        choices = ["RABIT", "LISA"], default = "LISA",
        help = "Method to predict driver regulators.")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")
    group_output.add_argument("--rseqc", dest = "rseqc", action = "store_true", 
        help = "Whether or not to run RSeQC. "
        "If set, the pipeline will include the RSeQC part and then takes a longer time. "
        "By default (not set), the pipeline will skip the RSeQC part.")

    return parser.parse_args()
    

def main():

    SCRIPT_PATH = os.path.dirname(__file__)

    myparser = CommandLineParser()
    directory = myparser.directory
    outpre = myparser.outprefix
    fastqdir = myparser.fastq_dir
    species = myparser.species
    platform = myparser.platform
    rseqc = myparser.rseqc
    method = myparser.method

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    if rseqc:
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scRNA_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

        readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_distr.png"%outpre)
        readqualplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_quality.png"%outpre)
        nvcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_NVC.png"%outpre)
        gcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_GCcontent.png"%outpre)
        genecovplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_genebody_cov.png"%outpre)
        countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)
        genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)
        geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)
        tdtitle = "Cluster-specific regulator identified by %s" %(method)
        tdcolname = "log10(%s score)" %(method)

        td_list = []
        for line in open(cluster_regulator_file,"r").readlines():
            if not line.startswith("Cluster"):
                items = line.strip().split("\t")
                items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
                items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
                td_list.append(items_str)
        td_str = "\n".join(td_list)

        #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
        report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "readdistr":readdistrplot_link,"readqual":readqualplot_link, "nvc":nvcplot_link, "gc":gcplot_link, "genecov":genecovplot_link, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtabletitle":tdtitle, "regtablecolname":tdcolname, "regtable":td_str}

    else:
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scRNA_noqc_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

        countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)
        genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)
        geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)
        tdtitle = "Cluster-specific regulator identified by %s" %(method)
        tdcolname = "log10(%s score)" %(method)

        td_list = []
        for line in open(cluster_regulator_file,"r").readlines():
            if not line.startswith("Cluster"):
                items = line.strip().split("\t")
                items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
                items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
                td_list.append(items_str)
        td_str = "\n".join(td_list)

        #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
        report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtabletitle":tdtitle, "regtablecolname":tdcolname, "regtable":td_str}


    report_html_instancefile = os.path.join(directory, outpre + "_scRNA_report.html")
    outf = open(report_html_instancefile,"w")
    outf.write(report_html_instance)
    outf.close()


if __name__ == "__main__":
    main()

