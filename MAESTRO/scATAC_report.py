#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 20:33:03 2019

@author: Dongqing Sun
"""

import os,sys

SCRIPT_PATH = os.path.dirname(__file__)
qcdir = sys.argv[1]
analysisdir = sys.argv[2]
outpre = sys.argv[3]

report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
report_html_temp = open(report_html_tempfile, "r").read()
fragplot_link = '''"%s/'''%qcdir + outpre + '''_scATAC_fragment_size.png"'''
mapplot_link = '''"%s/'''%qcdir + outpre + '''_scATAC_mapping_summary.png"'''
fripplot_link = '''"%s/'''%qcdir + outpre + '''_scATAC_cell_filtering.png"'''
peakcluster_link = '''"%s/'''%analysisdir + outpre + '''_cluster.png"'''
rpannotate_link = '''"%s/'''%analysisdir + outpre + '''_annotated.png"'''

report_html_instance = report_html_temp % {"fragment":fragplot_link, "map":mapplot_link, "frip":fripplot_link, "peakcluster":peakcluster_link, "rpannotate":rpannotate_link}
report_html_instancefile = "Result/" + outpre + "_scATAC_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()