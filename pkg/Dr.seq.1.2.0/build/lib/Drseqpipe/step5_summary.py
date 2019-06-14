#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string

# --------------------------
# custom package
# --------------------------

### tool function
from Drseqpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   CMD,
                                   createDIR,
                                   textformat,
                                   strlatexformat)
# --------------------------
# main 
# --------------------------
def step5_summary(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''
    # start
    # create section for 
    
    wlog('Step5: summary',logfile)
    wlog('copy results',logfile)
# Rscript analysis.r expmat outname coverGN highvarZ selectPCcutoff rdnumber maxKnum
    summarydir = conf_dict['General']['outputdirectory'] + 'summary/'
    createDIR(summarydir)
    os.chdir(summarydir)
    
    plot_folder = summarydir + "plots/"
    createDIR(plot_folder)
    os.chdir(plot_folder)
    ### collect results 
    for i in conf_dict['QCplots']:
        if os.path.isfile(conf_dict['QCplots'][i]):
            #realname
            cmd = 'cp %s .'%conf_dict['QCplots'][i]
            rwlog(cmd,logfile)

    result_folder = summarydir + "results/"
    createDIR(result_folder)
    os.chdir(result_folder)
    for i in conf_dict['results']:
        if os.path.isfile(conf_dict['results'][i]):
            cmd = 'cp %s .'%conf_dict['results'][i]
            rwlog(cmd,logfile)

    os.chdir(summarydir)

    wlog('generate qc documents',logfile)
    ### initiate 
    QCdoc = """\\documentclass[11pt,a4paper]{article}
\\usepackage{tabularx}
\\usepackage[english]{babel}
\\usepackage{array}
\\usepackage{graphicx}
\\usepackage{color}
\\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{QC and analysis reports for Drop-seq data : %s}

\\vspace{-1cm}
\\maketitle
\\tableofcontents
\\newpage
\\newpage
\\section{Data description}
\\begin{quotation}
Table 1 mainly describe the input file and mapping and analysis parameters.
\\end{quotation}
\\begin{table}[h]
\\caption{Data description}\\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(strlatexformat(conf_dict['General']['outname']))
    ### table1 prepare parameter
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        q30filter = "True"
    else:
        q30filter = "False"
    if int(conf_dict['Step2_ExpMat']['filterttsdistance']) == 1:
        filtertts = "True"
    else: 
        filtertts = "False"
    if int(conf_dict['Step2_ExpMat']['umidis1']) == 1:
        umidis1 = "True"
    else:
        umidis1 = "False"
    if int(conf_dict['Step3_QC']['remove_non_dup_cell']) == 1:
        rmnodup = "True"
    else:
        rmnodup = "False"
          
    QCdoc += """      
\\hline
parameter & value  \\\\
\\hline
output name & %s \\\\
\\hline
barcode file(file name only) & %s \\\\
\\hline
reads file(file name only) & %s \\\\
\\hline
reads file format & %s  \\\\
\\hline
cell barcode length &  %s \\\\
\\hline
UMI length & %s \\\\
\\hline
mapping software & %s \\\\
\\hline
Q30 filter mapped reads & %s \\\\
\\hline
remove reads away TTS & %s \\\\
\\hline
"""%(strlatexformat(conf_dict['General']['outname']),
     strlatexformat(conf_dict['General']['barcode_file'].split("/")[-1]),
     strlatexformat(conf_dict['General']['reads_file'].split("/")[-1]),
     conf_dict['General']['format'].upper(),
     str(conf_dict['General']['cell_barcode_length']),
     str(conf_dict['General']['umi_length']),
     conf_dict['Step1_Mapping']['mapping_software_main'],
     q30filter,
     filtertts
     )
    ### table1 part2
    if  filtertts == "True":
        QCdoc += """TTS distance (for remove) & %s bp \\\\
\\hline
"""%(str(conf_dict['Step2_ExpMat']['ttsdistance'])) 
    if  int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 1:
        QCdoc += """duplicate rate in each cell & UMI $+$ location \\\\"""
    elif int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 2:
        QCdoc += """duplicate rate in each cell & UMI only \\\\"""
    elif int(conf_dict['Step2_ExpMat']['duplicate_measure']) == 3:
        QCdoc += """duplicate rate in each cell & location only \\\\"""
    else:
        QCdoc += """duplicate rate in each cell & keep all reads \\\\"""
    if int(conf_dict['Step2_ExpMat']['duplicate_measure']) in [1,2]:
        QCdoc += """
\\hline
merge UMI ED = 1 & %s \\\\ 
\\hline"""%(umidis1)
    if  int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        QCdoc += """
select STAMPs & %s covered gene \\\\
\\hline"""%(str(conf_dict['Step3_QC']['covergncluster']))
    elif int(conf_dict['Step3_QC']['select_cell_measure']) == 2:
        QCdoc += """
select STAMPs & top %s UMI count \\\\
\\hline"""%(str(conf_dict['Step3_QC']['topumicellnumber']))
    QCdoc += """
remove low duplicate rate cell & %s \\\\ 
\\hline """%(rmnodup)
    if  rmnodup == "True":
        QCdoc += """
low duplicate rate cutoff & %s  \\\\
\\hline"""%(str(conf_dict['Step3_QC']['non_dup_cutoff']))
    QCdoc += """
z-score for highly variable gene & %s \\\\ 
\\hline 
cumulative variance for selecting PC & %s \\\\
\\hline """%(str(conf_dict['Step4_Analysis']['highvarz']),
     str(100*float(conf_dict['Step4_Analysis']['selectpccumvar']))+'\\%')
 
    if  int(conf_dict['Step4_Analysis']['clustering_method']) == 1:
        QCdoc += """
cluster method & k-means (Gap statistics, first stable) \\\\"""
    elif int(conf_dict['Step4_Analysis']['clustering_method']) == 2:
        QCdoc += """
cluster method & k-means (Gap statistics, maxSE) \\\\"""
    elif int(conf_dict['Step4_Analysis']['clustering_method']) == 3:
        QCdoc += """
cluster method & k-means (custom, k=%s) \\\\"""%(conf_dict['Step4_Analysis']['custom_k'])
    else:
        QCdoc += """
cluster method & DBScan (eps=%s) \\\\"""%(conf_dict['Step4_Analysis']['custom_d'])
    QCdoc += """
\\hline
\\end{tabularx}
\\end{table}
"""
    ### bulk QC
    QCdoc += """
\\newpage
\\newpage
\\section{Reads level QC}
In the reads level QC step we measured the quality of sequencing reads, including nucleotide quality and composition. In the reads level QC step and Bulk-cell level QC step we randomly sampled down total reads to 5 million and used a published package called ``RseQC" for reference.(Wang, L., Wang, S. and Li, W. (2012) )
\\newpage
\\newpage
\\subsection{Reads quality}
\\begin{quotation}
Reads quality is one of the basic reads level quality control methods. We plotted the distribution of a widely used Phred Quality Score at every position of sequence to measure the basic sequence quality of your data. Phred Quality Score was calculate by a python function $ord(Q) - 33$. Color in the heatmap represented frequency of this quality score observed at this position. Red represented higher frequency while blue was lower frequency. You may observe a decreasing of quality near the 3'end of sequence because of general degradation of quality over the duration of long runs. If the decreasing of quality influence the mappability (see ``Bulk-cell level QC") then the common remedy is to perform quality trimming where reads are truncated based on their average quality or you can trim serveal base pair near 3'end directly. If it doesn't help, you may consider your Drop-seq data poor quality. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads quality} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

\\newpage
\\newpage
\\subsection{Reads nucleotide composition}
\\begin{quotation}
We assess the nucleotide composition bias of a sample. The proportion of four different nucleotides was calculated at each position of reads. Theoretically four nucleotides had similar proportion at each position of reads. You may observe higher A/T count at 3'end of reads because of the 3'end polyA tail generated in sequencing cDNA libaray, otherwise the A/T count should be closer to C/G count. In any case, you should observe a stable pattern at least in the 3'end of reads. Spikes (un-stable pattern) which occur in the middle or tail of the reads indicate low sequence quality. You can trim serveral un-stable bases from the 3'end if low mappability (see ``Bulk-cell level QC") is also observed. If it doesn't help, you may consider your Drop-seq data poor quality. Note that t
he A/T vs G/C content can greatly vary from species to species. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads nucleotide composition} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

\\newpage
\\newpage
\\subsection{Reads GC content}
\\begin{quotation}
Distribution of GC content of each read. This module measures the general quality of the library. If the distribution looks different from a single bell (too sharp or too broad) then there may be a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species. If you observe sharp peak or broder peak and also observe low mappability (see ``Bulk-cell level QC"), you may consider your Drop-seq data poor quality.
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads GC content} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
"""%((conf_dict['QCplots']['read_qul'].split("/")[-1]),
     (conf_dict['QCplots']['read_nvc'].split("/")[-1]),
     (conf_dict['QCplots']['read_gc'].split("/")[-1])
    )

    QCdoc += """
\\newpage
\\newpage
\\section{Bulk-cell level QC}
In the bulk-cell level QC step we measured the performance of total Drop-seq reads. In this step we did't separate cell or remove ``empty" cell barcodes, just like treated the sample as bulk RNA-seq sample.
\\newpage
\\newpage
\\subsection{Reads alignment summary}
\\begin{quotation}
The following table shows mappability and distribution of total Drop-seq reads. It measures the general quality of data as a RNA-seq sample. Low mappability indicates poor sequence quality(see ``Reads level QC") or library quality(caused by contaminant). High duplicate rate (low total UMI percentage observed, e.g. $<$ 10\\%%) indicate insufficient RNA material and Overamplification. In summary, if the percentage of ``total UMI count" is less than 5\\%%, users may consider reconstruct your library(redo the experiment), but first you should make sure you already trim the adapter and map your reads to the corresponded species(genome version). Note that UMI number was calculated by removing duplicate reads (which have identical genomic location, cell barcode and UMI sequences). Mappable reads was after Q30 filtering if Q30 filter function was turned on.\\\\
** the percentage was calculated by dividing total reads number \\\\
*** the percentage was calculated by divding total UMI number
\\end{quotation}
\\begin{table}[h]
\\caption{Reads alignment summary}\\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|X| }
    
\\hline
genomic region(Category) &  reads number \\\\
\\hline
total reads & %s \\\\
\\hline
mappble reads &  %s (%s\\%%)* \\\\
\\hline
total UMI count & %s (%s\\%%)* \\\\
\\hline
CDS exon UMI count & %s (%s\\%%)** \\\\
\\hline
3'UTR UMI count & %s (%s\\%%)** \\\\
\\hline
5'UTR UMI count & %s (%s\\%%)** \\\\
\\hline
intron UMI count & %s (%s\\%%)** \\\\
\\hline
intergenic UMI count & %s (%s\\%%)** \\\\
\\hline

\\end{tabularx}
\\end{table}
"""%(textformat(str(conf_dict['Mapping_stat']['totalreads'])),
     textformat(str(conf_dict['Mapping_stat']['q30reads'])),
     str( round(100*conf_dict['Mapping_stat']['q30reads']*1.0/conf_dict['Mapping_stat']['totalreads'], 2)),
     textformat(str(conf_dict['Mapping_stat']['umi_gene'])),
     str( round(100*conf_dict['Mapping_stat']['umi_gene']*1.0/conf_dict['Mapping_stat']['totalreads'], 2)),
     textformat(str(conf_dict['Mapping_stat']['cdsN'])),
     str( round(100*conf_dict['Mapping_stat']['cdsN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['utr3N'])),
     str( round(100*conf_dict['Mapping_stat']['utr3N']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['utr5N'])),
     str( round(100*conf_dict['Mapping_stat']['utr5N']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['intronN'])),
     str( round(100*conf_dict['Mapping_stat']['intronN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)),
     textformat(str(conf_dict['Mapping_stat']['intergenicN'])),
     str( round(100*conf_dict['Mapping_stat']['intergenicN']*1.0/conf_dict['Mapping_stat']['umi_gene'], 2)))
     ### genebody coverage
    QCdoc += """
\\newpage
\\newpage
\\subsection{Gene body coverage}
\\begin{quotation}
Aggregate plot of reads coverage on all genes. This module measures the general quality of the Drop-seq data. Theoretically we observe a unimodal (single bell) distribution, but for Drop-seq sample an enrichment at 3'end is observed due to library preparation using oligo-dT primers. In any case you should observe a smooth distritbuion. If loss of reads or spike are observed in certain part of gene body (e.g. middle or 3'end of gene body), poor quality of your library was indicated. Especially when low mappability and high intron rate are also observed (see ``Reads alignment summary" section).
\\end{quotation}
\\begin{figure}[h]
        \\caption{Gene body coverage} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
"""%((conf_dict['QCplots']['gb_cover'].split("/")[-1]))

    QCdoc += """

\\newpage
\\newpage
\\section{Individual-cell level QC}
In this step we focused on the quality of individual cell and distinguishing cell barcodes from STAMPs (single-cell transcriptomes attached to microparticles)
\\newpage
\\newpage
\\subsection{Reads duplicate rate distribution}
\\begin{quotation}
Drop-seq technology has an innate advantage of detecting duplicate reads and amplification bias due to the barcode and UMI information. This module displays the distribution of duplicate rate in each cell barcode and helps to discard barcodes with low duplicate rate (which usually caused by empty cell barcodes and ambient RNA). We plot the distribution of duplicate rate in each cell barcode (though most of cell barcodes don't contain cells, they still have RNA) and observed a bimodal distribution of duplicate rate. We set an option for you to discard cell barcodes with low duplicate rate in following steps. The vertical line represented the cutoff (duplicate rate $>=$ 0.1) of discarding cell barcodes with low duplicate rate. You can adjust the cutoff and rerun Dr.seq if current cutoff didn't separate two peaks from the distribution clearly (usually happened with insufficient sequencing depth). If the distribution didn't show clear bimodal or you don't want to discard cell barcodes according to duplicate rate, you can set cutoff to 0 to keep all cell barcodes for following steps. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads dupliate rate distribution} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
"""%(conf_dict['QCplots']['duprate'].split("/")[-1])
    if int(conf_dict['Step3_QC']['select_cell_measure']) == 1:
        QCdoc += """
\\newpage
\\newpage
\\subsection{Reads duplicate rate vs. cumulative covered gene number}
\\begin{quotation}
Reads duplicate rate versus cumulative covered gene numbers. This module measures whether each of your individual cell was sequenced and clearly separated from empty cell barcodes. Cell barcodes are ranked by the number of covered genes. The duplicate rate (y-axis, left side) is plotted as a function of ranked cell barcode. Red curve represents the number of genes covered by top N cell barcodes (y-axis, right side). N is displayed by x-axis. Theoretically you observe a ``knee" on your cumulative curve (slope $=$ 1 on the curve) and the cutoff of your selected STAMPs (dash line) should be close to the ``knee". The cutoff can also be far away from the ``knee" in some cases because you input too many cells and have insufficient average sequencing depth, then you should adjust your cutoff (to the position you get enough STAMPs and sufficient reads count) and rerun Dr.seq. See the description of the paramter ``select cell measure" in the Manual.
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads duplicate rate vs. cumulative covered gene number} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

\\newpage
\\newpage
\\subsection{UMI vs. covered gene number}
\\begin{quotation}
Covered gene number is plotted as a function of the number of UMI (i.e. unique read). This module measures the quality of Drop-seq experiment and helps to distinguish STAMPs from empty cell barcodes. We observe a clearly different pattern for two groups of cell barcodes with different reads duplicate rate (blue dots versus red and purple dots). Purple dots represented the selected STAMPs for the cell-clustering step. By default we select STAMPs with 1000 gene covered after discarding low duplicate cell barcodes. You may get few STAMPs according to this cutoff if the average sequencing depth of your cells was too low or too many cells were inputed. In this case you can adjust your cutoff or tell Dr.seq to directly select cell barcodes with highest reads count (see the description of the parameter ``select cell measure"). Note that we use only STAMPs selected in this step for following analysis. The other cell barcodes are discarded. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{UMI v.s. covered gene number} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
"""%(conf_dict['QCplots']['cumumiduprate'].split("/")[-1],
     conf_dict['QCplots']['umicovergn'].split("/")[-1])
    else:
        QCdoc += """
\\newpage
\\newpage
\\subsection{Reads duplicate rate vs. cumulative covered gene number}
\\begin{quotation}
Reads duplicate rate versus cumulative covered gene numbers. This module measures whether each of your individual cell was sequenced and clearly separated from empty cell barcodes. Cell barcodes are ranked by the number of UMI count. The duplicate rate (y-axis, left side) is plotted as a function of ranked cell barcode. Red curve represents the number of genes covered by top N cell barcodes (y-axis, right side). N is displayed by x-axis. Theoretically you observe a ``knee" on your cumulative curve (slope $=$ 1 on the curve) and the cutoff of your selected STAMPs (dash line) should be close to the ``knee". The cutoff can also be far away from the ``knee" in some cases because you input too many cells and have insufficient average sequencing depth, then you should adjust your cutoff (to the position you get enough STAMPs and sufficient reads count) and rerun Dr.seq. See the description of the paramter ``select cell measure" in the Manual.
\\end{quotation}
\\begin{figure}[h]
        \\caption{Reads duplicate rate vs. cumulative covered gene number} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

\\newpage
\\newpage
\\subsection{UMI vs. covered gene number}
\\begin{quotation}
Covered gene number is plotted as a function of the number of UMI (i.e. unique read). This module measures the quality of Drop-seq experiment and helps to distinguish STAMPs from empty cell barcodes. We observe a clearly different pattern for two groups of cell barcodes with different reads duplicate rate (blue dots versus red and purple dots). Purple dots represented the selected STAMPs for the cell-clustering step. We select 1000 STAMPs with highest UMI count after discarding low duplicate cell barcodes. You may get few STAMPs according to this cutoff if the average sequencing depth of your cells was too low or too many cells were inputed. Note that we use only STAMPs selected in this step for following analysis. The other cell barcodes are discarded. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{UMI v.s. covered gene number} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
"""%(conf_dict['QCplots']['cumumiduprate'].split("/")[-1],
     conf_dict['QCplots']['umicovergn'].split("/")[-1])
     
    QCdoc += """
\\newpage
\\newpage
\\subsection{Covered gene number distribution}
\\begin{quotation}
Histogram of covered gene number of selected STAMPs. The module measures whether the selected STAMPs have sufficient reads coverage. By default Dr.seq selects cell barcodes with $>=$ 1000 genes covered as STAMPs. If you choose to select STAMPs with highest reads count (``select cell measure" $=$ 2), then you should check this figure to make sure the STAMPs you select have enough gene covered. If most of your STAMPs have low covered gene number (e.g. $<$ 100 gene covered), you can make your cutoff more stringent (e.g. select less cell barcodes with higher reads count) to make sure you get reliable STAMPs.
\\end{quotation}
\\begin{figure}[h]
        \\caption{Covered gene number} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

\\newpage
\\newpage
\\subsection{Intron rate distribution}
\\begin{quotation}
Intron rate is a effective method to measure the quality of a RNA-seq sample. We plot a histogram of intron rate of every STAMP barcodes to check whether reads from each STAMPs enriched in the exon region. High intron rate (e.g. $>=$ 30\\%%) indicates low quality of RNA in each STAMPs (caused by different problem, for example contaminant). You may consider your Drop-seq data low quality if most of selected STAMPs have high intron rate and low covered gene number (see ``Covered gene number distribution" section). Intron rate is defined as $\\frac{intron\\ reads\\ number}{intron + exon\\ reads\\ number}$ 
\\end{quotation}
\\begin{figure}[h]
        \\caption{Intron rate distribution} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
 
"""%(conf_dict['QCplots']['covergn'].split("/")[-1],
     conf_dict['QCplots']['intronrate'].split("/")[-1])
    
    if int(conf_dict['Step4_Analysis']['clustering_method']) in [3,4]:
        pass
    else:
        if int(conf_dict['Step4_Analysis']['clustering_method']) in [1,2]:
            selectM = 'first stable gap'
        else:
            selectM = 'maxSE'
        QCdoc += """
\\newpage
\\newpage
\\section{Cell-clustering level QC}
This step composed by k-means clustering based on t-SNE dimentional reduction result and Gap statistics to determine best k.
\\newpage
\\newpage
\\subsection{Gap statistics}
\\begin{quotation}
We conducted a k-means clustering based on t-SNE dimensional reduction output to measure sample's ability to be separated to different cell subtypes. Gap statistics was performed to determine the best k in k-means clustering. In general, decreasing pattern (usually k $<=$ 2) is observed for pure cell type or cell line data, while increasing pattern with bigger k should be observed for mix cell types (or cell subtypes) data. If the cluster number predicted from the Gap statistics is largely different to what you expect, it indicated that your cells are not well characterized and separated by the Drop-seq experiment (due to the contaminant or the low capture efficiency of Droplets). In this case, you may consider your Drop-seq data poor quality. Alternatively, you may would like to use the parameter ``custom k" to specify the cluster number.
\\end{quotation}
\\begin{figure}[h]
        \\caption{Gap statistics} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}

"""%(conf_dict['QCplots']['gapstat'].split("/")[-1])
    
    QCdoc += """
\\newpage
\\newpage
\\subsection{Clustering plot}
\\begin{quotation}
Scatter plot represented visualization of t-SNE dimensional reduction output of selected STAMP barcodes. STAMP barcodes are colored according to the clustering result and cluster numbers are printed in the center of each cluster. This figure is mainly for visualization and help you to know how your Drop-seq data look like. If you want to combine some small groups which are close to each other, you can use the cluster matrix (named ``cluster.txt") in the Dr.seq standard analysis output to conduct your own analysis.   
\\end{quotation}
\\begin{figure}[h]
        \\caption{Clustering plot} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
 
"""%(conf_dict['QCplots']['cluster'].split("/")[-1])
    if os.path.isfile(conf_dict['QCplots']['silhouette']):
        QCdoc += """
\\newpage
\\newpage
\\subsection{Silhouette of clustering}
\\begin{quotation}
Silhouette method is used to interprate and validate the consistency within clusters defined in previous steps. A poor Silhouette (e.g. average si $<$ 0.2 ) score indicate that Drop-seq experiments(if not properly done) may not separate well the subpopulations of cells. If most of your clusters have poor Silhouette score, it may indicate a poor quality of your Drop-seq experiments. 
\\end{quotation}
\\begin{figure}[h]
        \\caption{Silhouette score for clustered STAMPs} \\label{fig:profileunion}
        \\setlength{\\abovecaptionskip}{0pt}
        \\setlength{\\belowcaptionskip}{10pt}
        \\centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\\end{figure}
 
"""%(conf_dict['QCplots']['silhouette'].split("/")[-1])
    
#    QCdoc += """
#\\newpage
#\\newpage
#\subsection{STAMPs colored by total UMI count}
#\\begin{quotation}
#STAMPs was by the total number of UMI based on t-SNE visualization. 
#\end{quotation}
#\\begin{figure}[h]
#        \caption{STAMPs colored by total UMI count} \label{fig:profileunion}
#        \setlength{\\abovecaptionskip}{0pt}
#        \setlength{\\belowcaptionskip}{10pt}
#        \centering
#        {\includegraphics[width=0.8\\textwidth]{%s}}
#\end{figure}
 
#"""%(conf_dict['QCplots']['umicolor'].split("/")[-1])
   
#    QCdoc += """
#\\newpage
#\\newpage
#\subsection{STAMPs colored by intron rate}
#\\begin{quotation}
#STAMPs was by the intron rate based on t-SNE visualization. 
#\end{quotation}
#\\begin{figure}[h]
#        \caption{STAMPs colored by intron rate} \label{fig:profileunion}
#        \setlength{\\abovecaptionskip}{0pt}
#        \setlength{\\belowcaptionskip}{10pt}
#        \centering
#        {\includegraphics[width=0.8\\textwidth]{%s}}
#\end{figure}
# 
#"""%(conf_dict['QCplots']['itrcolor'].split("/")[-1])
      
    QCdoc += """
\\newpage
\\newpage
\\section{Output list}
\\begin{quotation}
All output files were described in the following table
\\end{quotation}
\\begin{table}[h]
\\caption{output list}\\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }
    
\\hline
description & filename \\\\
\\hline
expression matrix for selected STAMPs & %s  \\\\
"""%(strlatexformat(conf_dict['results']['expmatcc'].split("/")[-1]))
    if int(conf_dict['Step4_Analysis']['pctable']) == 1:
        QCdoc += """
\\hline
top2 components of PCA dimentional reduction result & %s \\\\         
"""%(strlatexformat(conf_dict['results']['pctable'].split("/")[-1]))
    if int(conf_dict['Step4_Analysis']['cortable']) == 1:
        QCdoc += """
\\hline
pairwise correlation matrix & %s \\\\
"""%(strlatexformat(conf_dict['results']['cortable'].split("/")[-1]))
    QCdoc += """
\\hline
All features of selected STAMPs & %s \\\\
\\hline
summary QC report & %s \\\\
\\hline

\\end{tabularx}
\\end{table} 
\\end{document} 
"""%(strlatexformat(conf_dict['results']['features'].split("/")[-1]),strlatexformat(conf_dict['General']['outname'])+"\\_summary.pdf")

    os.chdir(plot_folder)

    latexfile = conf_dict['General']['outname'] + '_summary.tex'
    outf = open(latexfile,'w')
    outf.write(QCdoc)
    outf.close()
    cmd = "pdflatex %s"%(latexfile)
    cmd2 = 'cp %s %s'%(conf_dict['General']['outname'] + '_summary.pdf',summarydir)
    if conf_dict['General']['latex'] == 1:
        rwlog(cmd,logfile)
        rwlog(cmd,logfile)
        rwlog(cmd2,logfile)
        for files in os.listdir(plot_folder):
            if os.path.isfile(files) and files[-12:-4] == "_summary":
                if not files[-4:] in ['.tex','.pdf',',png','.txt']:
                    cmd = "rm %s"%(files)
                    rwlog(cmd,logfile)
        wlog('pdflatex was detected in default PATH, generate summary report %s'%('summary/'+conf_dict['General']['outname'] + '_summary.pdf'),logfile)
    else:
        wlog('pdflatex was not detected in default PATH, generate summary report .tex file in summary/plots folder, you can move the whole summary/plots/ folder to the environment with pdflatex installed and run cmd in the plots/ folder: "pdflatex %s"'%(conf_dict['General']['outname'] + '_summary.tex'),logfile)
   
        
    if conf_dict['clean']:
        wlog('--clean pararmeter was turned on, remove internal files with large size',logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_symbol.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_cds.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_3utr.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_5utr.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_on_TTSdis.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_combined.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'expmatrix/' + conf_dict['General']['outname']+'_barcode_reform.txt'),logfile)

    wlog('Step5 summary DONE, check %s for final outputs'%(summarydir),logfile)


    return conf_dict









