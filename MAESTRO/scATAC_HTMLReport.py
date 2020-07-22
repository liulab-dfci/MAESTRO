# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-28 03:15:37
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-07-22 15:59:35


import os,sys
import snakemake.report

import argparse as ap

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Generate scATAC result report. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--platform", dest = "platform", default = "10x-genomics", 
        choices = ["10x-genomics", "sci-ATAC-seq", "microfluidic"], 
        help = "Platform of single cell ATAC-seq. DEFAULT: 10x-genomics.")
    group_input.add_argument("--input-format", dest = "input_format", default = "fastq", 
        choices = ["fastq", "bam", "fragments"], 
        help = "The format of input files. DEFAULT: fastq.")
    group_input.add_argument("--input-path", dest = "input_path", type = str, default = "",  
        help = "Directory where input files are stored")
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
    input_path = myparser.input_path
    species = myparser.species
    platform = myparser.platform
    input_format = myparser.input_format

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    fragplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_fragment_size.png"%outpre)
    # mapplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_mapping_summary.png"%outpre)[0]
    fripplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_cell_filtering.png"%outpre)
    peakcluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)
    if os.path.exists("Result/Analysis/%s_annotated.png"%outpre):
        rpannodisplay = "inline"
        rpannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)
    else:
        rpannodisplay = "none"
        rpannotate_link = ""
    if os.path.exists("Result/Analysis/%s_CistromeTop_annotated.png"%outpre):
        caannodisplay = "inline"
        caannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_CistromeTop_annotated.png"%outpre)
    else:
        caannodisplay = "none"
        caannotate_link = ""

    if os.path.exists("Result/Analysis/%s_MS4A1_genetrack.png"%outpre):
        ms4a1display = "inline"
        ms4a1track_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_MS4A1_genetrack.png"%outpre)
    else:
        ms4a1display = "none"
        ms4a1track_link = ""

    if os.path.exists("Result/Analysis/%s_CD3D_genetrack.png"%outpre):
        cd3ddisplay = "inline"
        cd3dtrack_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_CD3D_genetrack.png"%outpre)
    else:
        cd3ddisplay = "none"
        cd3dtrack_link = ""

    td_list = []
    for line in open(cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            td_list.append(items_str)
    td_str = "\n".join(td_list)

    if input_format != "fragments":
        if input_format == "fastq":
            inputformat = "FASTQ Path"
        else:
            inputformat = "BAM Path"
        readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_read_distr.png"%outpre)
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        report_html_instance = report_html_temp % {"outprefix":outpre, "inputformat":inputformat, "inputpath":input_path, "species":species, 
        "platform":platform, "readdistr":readdistrplot_link, "fragment":fragplot_link, "frip":fripplot_link, 
        "peakcluster":peakcluster_link, "rpannodisplay": rpannodisplay, "rpannotate":rpannotate_link, "regtable":td_str,
        "caannodisplay": caannodisplay, "caannotate": caannotate_link, 'cd3ddisplay': cd3ddisplay, 'ms4a1display': ms4a1display,
        'cd3dtrack': cd3dtrack_link, "ms4a1track": ms4a1track_link}
    else:
        inputformat = "Fragment Path"
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_nomappability_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        report_html_instance = report_html_temp % {"outprefix":outpre, "inputformat":inputformat, "inputpath":input_path, "species":species, 
        "platform":platform, "fragment":fragplot_link, "frip":fripplot_link, 
        "peakcluster":peakcluster_link, "rpannodisplay": rpannodisplay, "rpannotate":rpannotate_link, "regtable":td_str,
        "caannodisplay": caannodisplay, "caannotate": caannotate_link, 'cd3ddisplay': cd3ddisplay, 'ms4a1display': ms4a1display,
        'cd3dtrack': cd3dtrack_link, "ms4a1track": ms4a1track_link}

    report_html_instancefile = os.path.join(directory, outpre + "_scATAC_report.html")
    outf = open(report_html_instancefile,"w")
    outf.write(report_html_instance)
    outf.close()


if __name__ == "__main__":
    main()

