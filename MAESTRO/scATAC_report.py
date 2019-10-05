#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:33:03 2019

@author: Dongqing Sun
"""

import os,sys
import snakemake.report

SCRIPT_PATH = os.path.dirname(__file__)
outpre = sys.argv[1]
fastqdir = sys.argv[2]
species = sys.argv[3]
platform = sys.argv[4]

report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
report_html_temp = open(report_html_tempfile, "r").read()

# readdistrplot_link = '''"Plot/%s_scATAC_read_distr.png"'''%outpre
# fragplot_link = '''"Plot/%s_scATAC_fragment_size.png"'''%outpre
# mapplot_link = '''"Plot/%s_scATAC_mapping_summary.png"'''%outpre
# fripplot_link = '''"Plot/%s_scATAC_cell_filtering.png"'''%outpre
# peakcluster_link = '''"Plot/%s_cluster.png"'''%outpre
# rpannotate_link = '''"Plot/%s_annotated.png"'''%outpre
bulkqc_file = "Result/QC/%s_bam_stat.txt"%outpre
cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre



fragplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_fragment_size.png"%outpre)[0]
mapplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_mapping_summary.png"%outpre)[0]
fripplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_cell_filtering.png"%outpre)[0]
peakcluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)[0]
rpannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]
readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scATAC_read_distr.png"%outpre)[0]

line_id = 0
total,mapped,duplicate,mito,uniq,promoters = 0,0,0,0,0,0
for line in open(bulkqc_file).readlines():
    line = line.strip().split(' ')
    line_id += 1
    if line_id == 1:
        total = int(line[0])
    if line_id == 5:
        mapped = int(line[0])
    if line_id == 4:
        duplicate = int(line[0])
    if line_id == 14:
        mito = int(line[0])
    if line_id == 15:
        uniq = int(line[0])
    if line_id == 16:
        promoters = int(line[0])
stat_list = [total,mapped,duplicate,mito,uniq,promoters]

stat_list[1] = str(stat_list[1]) + "(%.2f" %(float(stat_list[1])/float(stat_list[0])*100) + "%)"
stat_list[2] = str(stat_list[2]) + "(%.2f" %(float(stat_list[2])/float(stat_list[0])*100) + "%)"
stat_list[3] = str(stat_list[3]) + "(%.2f" %(float(stat_list[3])/float(stat_list[0])*100) + "%)"
stat_list[4] = str(stat_list[4]) + "(%.2f" %(float(stat_list[4])/float(stat_list[0])*100) + "%)"
stat_list[5] = str(stat_list[5]) + "(%.2f" %(float(stat_list[5])/float(stat_list[0])*100) + "%)"

td_list = []
for line in open(cluster_regulator_file,"r").readlines():
    if not line.startswith("Cluster"):
        items = line.strip().split("\t")
        items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
        items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
        td_list.append(items_str)
td_str = "\n".join(td_list)

report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species, "platform":platform, "readdistr":readdistrplot_link, "fragment":fragplot_link, "frip":fripplot_link, "peakcluster":peakcluster_link, "rpannotate":rpannotate_link, "regtable":td_str}

report_html_instancefile = "Result/" + outpre + "_scATAC_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()