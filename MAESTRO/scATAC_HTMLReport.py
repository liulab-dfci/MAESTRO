# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-28 03:15:37
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-02-29 04:07:31


import os,sys
import snakemake.report

import argparse as ap

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Generate scATAC result report. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--platform", dest = "platform", default = "10x-genomics", 
        choices = ["10x-genomics", "sci-ATAC-seq", "microfluidic"], 
        help = "Platform of single cell ATAC-seq. DEFAULT: 10x-genomics.")
    group_input.add_argument("--fastq-dir", dest = "fastq_dir", type = str, default = "",  
        help = "Directory where fastq files are stored")
    group_input.add_argument("--species", dest = "species", default = "GRCh38", 
        choices = ["GRCh38", "GRCm38"], type = str, 
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics", 
        help = "Prefix of output files. DEFAULT: MAESTRO.")
    
    return parser.parse_args()
    

def main():

    SCRIPT_PATH = os.path.dirname(__file__)

    myparser = CommandLineParser()
    directory = myparser.directory
    outpre = myparser.outprefix
    fastqdir = myparser.fastq_dir
    species = myparser.species
    platform = myparser.platform

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
    report_html_temp = open(report_html_tempfile, "r").read()

    cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    fragplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_fragment_size.png"%outpre)
    # mapplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_mapping_summary.png"%outpre)[0]
    fripplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_cell_filtering.png"%outpre)
    peakcluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)
    rpannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)
    readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_read_distr.png"%outpre)


    td_list = []
    for line in open(cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            td_list.append(items_str)
    td_str = "\n".join(td_list)

    report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species, "platform":platform, "readdistr":readdistrplot_link, "fragment":fragplot_link, "frip":fripplot_link, "peakcluster":peakcluster_link, "rpannotate":rpannotate_link, "regtable":td_str}

    report_html_instancefile = os.path.join(directory, outpre + "_scATAC_report.html")
    outf = open(report_html_instancefile,"w")
    outf.write(report_html_instance)
    outf.close()


if __name__ == "__main__":
    main()

