#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:33:03 2019

@author: Dongqing Sun
"""

import os,sys

SCRIPT_PATH = os.path.dirname(__file__)
outpre = sys.argv[1]

report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
report_html_temp = open(report_html_tempfile, "r").read()
fragplot_link = '''"QC/''' + outpre + '''_frag.png"'''
mapplot_link = '''"QC/''' + outpre + '''_map.png"'''
fripplot_link = '''"QC/''' + outpre + '''_frip.png"'''
peakcluster_link = '''"Analysis/''' + outpre + '''_UMAP_cluster.png"'''
rpannotate_link = '''"Analysis/UMAP_assignIdent_''' + outpre + '''_RPannotated.png"'''

report_html_instance = report_html_temp % {"fragment":fragplot_link, "map":mapplot_link, "frip":fripplot_link, "peakcluster":peakcluster_link, "rpannotate":rpannotate_link}
report_html_instancefile = "Result/" + outpre + "_scATAC_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()