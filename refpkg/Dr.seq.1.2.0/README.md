Time-stamp: <2015-12-07 Shengen Hu> Tarelahu@gmail.com

webpage: 
         
http://www.tongji.edu.cn/~zhanglab/drseq 

Introduction: Drop-seq has recently emerged as a powerful technology to analyze gene expression from thousands of individual cells simultaneously. Currently, Drop-seq technology requires refinement, and quality control (QC) steps are critical for such data analysis. There is a strong need for a convenient and comprehensive approach to obtain dedicated QC and to determine the relationships between cells for ultra-high-dimensional data sets. We developed Dr.seq, a QC and analysis pipeline for Drop-seq data. By applying this pipe-line, Dr.seq provides four groups of QC measurements for given Drop-seq data, including reads level, bulk-cell level, individual-cell level and cell-clustering level QC. We assessed Dr.seq on simulated and published Drop-seq data. Both assessments exhibit reliable results. Overall, Dr.seq is a compre-hensive QC and analysis pipeline designed for Drop-seq data that is easily extended to other droplet-based data types.
Main propose: Dr.seq can generate detailed quality control report and analysis result of a Drop-seq sample, with 2 FASTQ files from a single sample inputted (reads file and barcode file)

Check documents in README/ folder for quick start and usage.

Quick start: Get you start and familiar with Dr.seq within 3 steps,
    1.Install,
    2.Download annoatation,
    3.Run

Manual: Full manual including
    1.Installation of Dr.seq,
    2.Preparation of Dr.seq related data and package,
    3.Usage of 2 mode of Dr.seq,
    4.Description of all changeable parameter

FAQ: frequently asked question (Email me if you have additional question, I'll add to the FAQ list.)

Changelog v1.2: set default of measure duplicate as 2 (only consider UMI)
