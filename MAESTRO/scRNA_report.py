#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 15:59:29 2019

@author: Dongqing Sun
"""

import os,sys
import snakemake.report

SCRIPT_PATH = os.path.dirname(__file__)
outpre = sys.argv[1]
fastqdir = sys.argv[2]
species = sys.argv[3]
platform = sys.argv[4]
rseqc = sys.argv[5]


if rseqc == "True":
    report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scRNA_template.html")
    report_html_temp = open(report_html_tempfile, "r").read()

    cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_distr.png"%outpre)[0]
    readqualplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_quality.png"%outpre)[0]
    nvcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_NVC.png"%outpre)[0]
    gcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_GCcontent.png"%outpre)[0]
    genecovplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_genebody_cov.png"%outpre)[0]
    countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)[0]
    genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)[0]
    geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]

    td_list = []
    for line in open(cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            td_list.append(items_str)
    td_str = "\n".join(td_list)

    #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
    report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "readdistr":readdistrplot_link,"readqual":readqualplot_link, "nvc":nvcplot_link, "gc":gcplot_link, "genecov":genecovplot_link, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtable":td_str}

elif platform == "10x-genomics":
    report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scRNA_10x_noqc_template.html")
    report_html_temp = open(report_html_tempfile, "r").read()

    cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_distr.png"%outpre)[0]
    countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)[0]
    genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)[0]
    geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]

    td_list = []
    for line in open(cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            td_list.append(items_str)
    td_str = "\n".join(td_list)

    #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
    report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "readdistr":readdistrplot_link, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtable":td_str}


else:
    report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scRNA_noqc_template.html")
    report_html_temp = open(report_html_tempfile, "r").read()

    cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre

    countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)[0]
    genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)[0]
    geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]

    td_list = []
    for line in open(cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            items = line.strip().split("\t")
            items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
            items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
            td_list.append(items_str)
    td_str = "\n".join(td_list)

    #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
    report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtable":td_str}


# readdistrplot_link = '''"Plot/%s_scRNA_read_distr.png"'''%outpre
# readqualplot_link = '''"Plot/%s_scRNA_read_quality.png"'''%outpre
# nvcplot_link = '''"Plot/%s_scRNA_NVC.png"'''%outpre
# gcplot_link = '''"Plot/%s_scRNA_GCcontent.png"'''%outpre
# genecovplot_link = '''"Plot/%s_scRNA_genebody_cov.png"'''%outpre
# countgeneplot_link = '''"Plot/%s_scRNA_cell_filtering.png"'''%outpre
# genecluster_link = '''"Plot/%s_cluster.png"'''%outpre
# geneannotate_link = '''"Plot/%s_annotated.png"'''%outpre

# bamstat_file = "Result/QC/%s_bam_stat.txt"%outpre
# readdistr_file = "Result/QC/%s_read_distribution.txt"%outpre
# cluster_regulator_file = "Result/Analysis/%s.PredictedTFTop10.txt"%outpre


# readdistrplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_distr.png"%outpre)[0]
# readqualplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_read_quality.png"%outpre)[0]
# nvcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_NVC.png"%outpre)[0]
# gcplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_GCcontent.png"%outpre)[0]
# genecovplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_genebody_cov.png"%outpre)[0]
# countgeneplot_link = snakemake.report.data_uri_from_file("Result/QC/%s_scRNA_cell_filtering.png"%outpre)[0]
# genecluster_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_cluster.png"%outpre)[0]
# geneannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/%s_annotated.png"%outpre)[0]

# stat_list = []
# for line in open(bamstat_file, "r").readlines():
#     if line.startswith("Total records:"):
#         total_reads = line.strip().split(":")[1].strip()
#         stat_list.append(int(total_reads))
#     if line.startswith("Optical/PCR duplicate:"):
#         dup_reads = line.strip().split(":")[1].strip()
#         stat_list.append(int(dup_reads))
#         break
#     else:
#         pass

# exon_tags = 0
# for line in open(readdistr_file, "r").readlines():
#     if line.startswith("Total Reads"):
#         mapped_reads = line.strip().split(" ")[-1]
#         stat_list.append(int(mapped_reads))
#         continue
#     if line.startswith("Total Tags"):
#         mapped_tags = line.strip().split(" ")[-1]
#         stat_list.append(int(mapped_tags))
#         continue
#     if line.startswith("CDS_Exons"):
#         cds_tags = line.strip().split()[2]
#         exon_tags += int(cds_tags)
#         continue
#     if line.startswith("5'UTR_Exons"):
#         utr5_tags = line.strip().split()[2]
#         exon_tags += int(utr5_tags)
#         continue
#     if line.startswith("3'UTR_Exons"):
#         utr3_tags = line.strip().split()[2]
#         exon_tags += int(utr3_tags)
#         stat_list.append(exon_tags)
#         continue
#     if line.startswith("Introns"):
#         intron_tags = line.strip().split()[2]
#         stat_list.append(int(intron_tags))
#         break

# # total reads, duplicate reads, mapped reads, mapped tags, exon tags, intron tags
# stat_list[1] = str(stat_list[1]) + "(%.2f" %(float(stat_list[1])/float(stat_list[0])*100) + "%)"
# stat_list[2] = str(stat_list[2]) + "(%.2f" %(float(stat_list[2])/float(stat_list[0])*100) + "%)"
# stat_list[4] = str(stat_list[4]) + "(%.2f" %(float(stat_list[4])/float(stat_list[3])*100) + "%)"
# stat_list[5] = str(stat_list[5]) + "(%.2f" %(float(stat_list[5])/float(stat_list[3])*100) + "%)"

# td_list = []
# for line in open(cluster_regulator_file,"r").readlines():
#     if not line.startswith("Cluster"):
#         items = line.strip().split("\t")
#         items_str_list = ["                                                            <td>" + i + "</td>" for i in items]
#         items_str = "                                                        <tr>\n" + "\n".join(items_str_list) + "\n                                                        </tr>"
#         td_list.append(items_str)
# td_str = "\n".join(td_list)

# #"totalreads":stat_list[0],"dupreads":stat_list[1],"mapreads":stat_list[2],"maptags":stat_list[3],"exontags":stat_list[4],"introntags":stat_list[5],
# report_html_instance = report_html_temp % {"outprefix":outpre, "fastqdir":fastqdir, "species":species,"platform":platform, "readdistr":readdistrplot_link,"readqual":readqualplot_link, "nvc":nvcplot_link, "gc":gcplot_link, "genecov":genecovplot_link, "countgene":countgeneplot_link, "genecluster":genecluster_link, "geneannotate":geneannotate_link, "regtable":td_str}

report_html_instancefile = "Result/" + outpre + "_scRNA_report.html"
outf = open(report_html_instancefile,"w")
outf.write(report_html_instance)
outf.close()
