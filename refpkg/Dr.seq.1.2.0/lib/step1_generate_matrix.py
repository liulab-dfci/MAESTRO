#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import time
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
                                   sample_down_transform_sam,
                                   transform_refgene,
                                   reform_barcode_fastq,
                                   combine_reads,
                                   generate_matrix
                                   )
# --------------------------
# main 
# --------------------------

def step1_generate_matrix(conf_dict,logfile):
    '''
    generate expression matrix file 
    main data processing step, including mapping, generate expression matrix and QC matrix which is used in next step
    for fastq format : 
        STAR/bowtie2 mapping
        q30 filter, 
    for sam format:
        q30 filter     
    ''' 
    wlog("Step1: alignment",logfile)
    t= time.time()
    ### create mapping dir 
    mapping_dir = conf_dict['General']['outputdirectory'] + 'mapping/'
    createDIR(mapping_dir)
    ### check reads file format , start mapping step if format is fastq
    if conf_dict['General']['format'] == 'sam':
        wlog('reads file format is sam, skip mapping step',logfile)
        conf_dict['General']['sam'] = conf_dict['General']['reads_file']
    else:
        wlog('Now start mapping in %s , all mapping result will be here'%(mapping_dir),logfile)
        os.chdir(mapping_dir)
        ## choose mapping tool from STAR and bowtie2 according to config file
        if conf_dict['Step1_Mapping']['mapping_software_main'] == "STAR":
            wlog('user choose STAR as alignment software',logfile)
            if sp('which STAR')[0].strip() == "":
                ewlog('STAR is not detected in default PATH, make sure you installed STAR and export it into default PATH',logfile)
            mapping_cmd = 'STAR --genomeDir %s --readFilesIn %s --runThreadN %s'%(conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['reads_file'],conf_dict['Step1_Mapping']['mapping_p'])
            mapping_cmd2 = 'mv Aligned.out.sam %s.sam'%(conf_dict['General']['outname'])
            rwlog(mapping_cmd,logfile)
            rwlog(mapping_cmd2,logfile)
            
        elif conf_dict['Step1_Mapping']['mapping_software_main'] == "bowtie2":
            wlog('user choose bowtie2 as alignment software',logfile)
            if sp('which bowtie2')[0].strip() == "":
                ewlog('bowtie2 is not detected in default PATH, make sure you installed bowtie2 and export it into default PATH',logfile)
            mapping_cmd = 'bowtie2 -p %s -x %s -U %s -S %s.sam   2>&1 >>/dev/null |tee -a %s.bowtieout'%(conf_dict['Step1_Mapping']['mapping_p'],conf_dict['Step1_Mapping']['mapindex'],conf_dict['General']['reads_file'],conf_dict['General']['outname'],conf_dict['General']['outname'])
            rwlog(mapping_cmd,logfile)
            
        else:
            ewlog("alignment tools can only be STAR and bowtie2",logfile)

        conf_dict['General']['sam'] = mapping_dir + conf_dict['General']['outname'] + '.sam'
    ### transform to bed file, awk helps to conduct q30 filtering
    wlog("transfer sam file to aligned bed file with own script",logfile)
    conf_dict['General']['bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed' 
    conf_dict['General']['sampledownsam'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.sam' 
    conf_dict['General']['sampledownbed'] = mapping_dir + conf_dict['General']['outname'] + '_sampledown.bed' 
    if int(conf_dict['Step1_Mapping']['q30filter']) == 1:
        wlog("q30 filter is turned on",logfile)
    else:
        wlog("q30 filter is turned off",logfile)
    ### use own script to transform sam to bed, and random sampling 5M mappable reads
    sample_down_transform_sam(conf_dict['General']['sam'],conf_dict['General']['bed'],conf_dict['General']['sampledownsam'],conf_dict['General']['sampledownbed'],5000000,int(conf_dict['Step1_Mapping']['q30filter']))
#        q30cmd = """samtools view -q 30 -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr" && $5 > 30) {if ($2 == 16) print $3,$4-1,$4,$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        rwlog(q30cmd,logfile,conf_dict['General']['dryrun'])
#        q30cmd = """samtools view -XS %s | awk '{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if (substr($2,1,1) == "r") print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4-1+length($11),$1,255,"-";else print $3,$4-1,$4-1+length($11),$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        q30cmd = """awk '/^[^@]/{FS="\t";OFS="\t";if (substr($3,1,3) == "chr") {if ($2 == 16) print $3,$4-1,$4+length($11),$1,255,"-";else print $3,$4-1,$4,$1,255,"+";}}' %s > %s"""%(conf_dict['General']['sam'],conf_dict['General']['bed'])
#        rwlog(q30cmd,logfile,conf_dict['General']['dryrun'])
    if not os.path.isfile(conf_dict['General']['bed']) or os.path.getsize(conf_dict['General']['bed']) == 0:
        ewlog('Alignment step / q30 filtering step failed, check your alignment parameter and samfile',logfile)
    s1time = time.time() -t
    wlog("time for alignment: %s"%(s1time),logfile)
    wlog("Step1: alignment DONE",logfile)

    ### create annotation dir and generate related annotation file
    t = time.time() 
    wlog("Step2: transform expression matrix",logfile)
    wlog('generate related annotation file with own script',logfile)
    annotation_dir = conf_dict['General']['outputdirectory'] + 'annotation/'
    createDIR(annotation_dir)
    os.chdir(annotation_dir)    
    transform_refgene(conf_dict['General']['gene_annotation'],conf_dict['Step2_ExpMat']['ttsdistance'],conf_dict['General']['outname'])

    ### create expression matrix dir and generate matrix
    wlog('generate expression matrix and individual cell qc matrix with own script',logfile)    
    expdir = conf_dict['General']['outputdirectory'] + 'expmatrix/'
    createDIR(expdir)
    os.chdir(expdir)
    
    ### use bedtools(intersect function) to assign exon/intron/intergenic/overlapping gene  information to all reads
    ### sort according to name
    wlog('add gene annotation on aligned bed file',logfile)
    cmd1 = "bedtools intersect -a %s -b %s  -wo   | sort -k 4,4 - >  %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_symbol.bed',conf_dict['General']['outname']+'_on_symbol.bed')
    cmd2 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_cds.bed',conf_dict['General']['outname']+'_on_cds.bed')
    cmd3 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_3utr.bed',conf_dict['General']['outname']+'_on_3utr.bed')
    cmd4 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_5utr.bed',conf_dict['General']['outname']+'_on_5utr.bed')
    cmd5 = "bedtools intersect -a %s -b %s -c | sort -k 4,4 - > %s"%(conf_dict['General']['bed'],annotation_dir+conf_dict['General']['outname']+'_gene_anno_TTSdis.bed',conf_dict['General']['outname']+'_on_TTSdis.bed')
    rwlog(cmd1,logfile)
    rwlog(cmd2,logfile)
    rwlog(cmd3,logfile)
    rwlog(cmd4,logfile)
    rwlog(cmd5,logfile)

    ### transform barcode fastq to 3column txt file [name,cell_barcode,umi]
    if conf_dict['General']['format1'] == 'txt':
        wlog('barcode files is reformed txt format, skip reform step',logfile)
        conf_dict['General']['barcode_reform'] = conf_dict['General']['barcode_file']
    else:
        wlog('reform barcode files with own script',logfile)
        conf_dict['General']['barcode_reform'] = expdir + conf_dict['General']['outname'] + '_barcode_reform.txt'
        reform_barcode_fastq(conf_dict['General']['barcode_file'],conf_dict['General']['barcode_reform'],conf_dict['General']['cell_barcode_length'],conf_dict['General']['umi_length'])
    ### sort according name
    cmdsort = 'sort -k 1,1 %s > %s'%(conf_dict['General']['barcode_reform'],expdir + conf_dict['General']['outname'] + '_barcode_reform_sort.txt')
    rwlog(cmdsort,logfile)
    conf_dict['General']['barcode_reform'] = expdir + conf_dict['General']['outname'] + '_barcode_reform_sort.txt'
    
    ### combine gene annotation, reads, barcode together
    wlog('combine annotation and barcode on reads with own script',logfile)
    combine_reads(conf_dict['General']['barcode_reform'],conf_dict['General']['outname']+'_on_cds.bed',conf_dict['General']['outname']+'_on_3utr.bed',conf_dict['General']['outname']+'_on_5utr.bed',conf_dict['General']['outname']+'_on_symbol.bed',conf_dict['General']['outname']+'_on_TTSdis.bed',conf_dict['General']['outname']+ '_combined.bed',conf_dict['Step2_ExpMat']['duplicate_measure'])   
     
    ### sort combined file by umi+loci, for following duplicate detection
    cmd6 = "sort -k 7,7 -k 5,5 %s > %s"%(conf_dict['General']['outname']+ '_combined.bed',conf_dict['General']['outname']+ '_combined_sort.bed')
    rwlog(cmd6,logfile)
    
    ### generate expression and QC matrix based on combined file
    wlog('generate expression matrix and QC matrix with own script',logfile)
    ### qcmatfull contains all cell_barcodes, while qcmat,expmat only contain cell_barcodes >= covergncutoff(100, default)
    conf_dict['Step2_ExpMat']['qcmatfull'] = expdir + conf_dict['General']['outname'] + "_qcmatfull.txt"    
    conf_dict['Step2_ExpMat']['qcmat'] = expdir + conf_dict['General']['outname'] + "_qcmat.txt"
    conf_dict['Step2_ExpMat']['expmat'] = expdir + conf_dict['General']['outname'] + "_expmat.txt"
    
    generate_matrix(conf_dict['General']['gene_annotation'],conf_dict['General']['outname']+ '_combined_sort.bed',conf_dict['Step2_ExpMat']['filterttsdistance'],conf_dict['Step2_ExpMat']['qcmatfull'],conf_dict['Step2_ExpMat']['qcmat'],conf_dict['Step2_ExpMat']['expmat'],conf_dict['Step2_ExpMat']['covergncutoff'],conf_dict['Step2_ExpMat']['umidis1'])    
        
    wlog("Step2 transform expression matrix DONE",logfile)
    s2time = time.time() -t
    wlog("time for transform expmat: %s"%(s2time),logfile)
    conf_dict['results'] = {}
    #conf_dict['results']['expmat'] = conf_dict['Step2_ExpMat']['expmat']
    #conf_dict['results']['qcmat'] = conf_dict['Step2_ExpMat']['qcmat']
    
    return conf_dict
