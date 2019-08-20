#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:33:03 2019

@author: Dongqing Sun
"""

import os,sys

SCRIPT_PATH = os.path.dirname(__file__)
outpre = sys.argv[1]
fastqdir = sys.argv[2]
species = sys.argv[3]
platform = sys.argv[4]

report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
report_html_temp = open(report_html_tempfile, "r").read()

fragplot_link = '''"Plot/''' + outpre + '''_scATAC_fragment_size.png"'''
mapplot_link = '''"Plot/''' + outpre + '''_scATAC_mapping_summary.png"'''
fripplot_link = '''"Plot/''' + outpre + '''_scATAC_cell_filtering.png"'''
peakcluster_link = '''"Plot/''' + outpre + '''_cluster.png"'''
rpannotate_link = '''"Plot/''' + outpre + '''_annotated.png"'''
cluster_regulator_file = "Result/Analysis/" + outpre + ".PredictedTFTop10.txt"

td_list = []
for line in open(cluster_regulator_file,"r").readlines():
    if not line.startswith("Cluster"):
        items = line.strip().split("\t")
        items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
        items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
        td_list.append(items_str)
td_str = "\n".join(td_list)

report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species, "platform":platform, "fragment":fragplot_link, "map":mapplot_link, "frip":fripplot_link, "peakcluster":peakcluster_link, "rpannotate":rpannotate_link, "regtable":td_str}

report_html_instancefile = "Result/Summary/" + outpre + "_scATAC_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()