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
rnaobj = sys.argv[2]
atacobj = sys.argv[3]

report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "Integration_template.html")
report_html_temp = open(report_html_tempfile, "r").read()

# allsourceplot_link = '''"Plot/''' + outpre + '''_source.png"'''
# rnaplot_link = '''"Plot/''' + outpre + '''_RNAonly.png"'''
# atacplot_link = '''"Plot/''' + outpre + '''_ATAConly.png"'''
# annoplot_link = '''"Plot/''' + outpre + '''_annotated.png"'''

allsourceplot_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_source.png"%outpre)[0]
rnaplot_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_RNAonly.png"%outpre)[0]
atacplot_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_ATAConly.png"%outpre)[0]
annoplot_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]

report_html_instance = report_html_temp % {"outprefix":outpre, "scrnaobject":rnaobj, "scatacaobject":atacobj, "allsource":allsourceplot_link, "rnacluster":rnaplot_link, "ataccluster":atacplot_link, "allannotate":annoplot_link}

report_html_instancefile = "Result/" + outpre + "_integrate_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()