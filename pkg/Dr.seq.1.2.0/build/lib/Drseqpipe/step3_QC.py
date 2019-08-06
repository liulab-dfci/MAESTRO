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
                                   readsqc,
                                   GBcover)
# --------------------------
# main 
# --------------------------
def step3_QC(conf_dict,logfile):
    '''
    start RseQC
    mapping stat
    single cell level QC
    '''
    # start
    # create section for 
    
    wlog('Step3: bulk and individual cell QC',logfile)
    ### preparing mapping state dict
    wlog('calculate mapping state',logfile)
    conf_dict['Mapping_stat'] = {}
    conf_dict['Mapping_stat']['umi_gene'] = 0
    conf_dict['Mapping_stat']['cdsN'] = 0
    conf_dict['Mapping_stat']['utr3N'] = 0
    conf_dict['Mapping_stat']['utr5N'] = 0
    conf_dict['Mapping_stat']['intronN'] = 0
    conf_dict['Mapping_stat']['intergenicN'] = 0

    ### calculate mapping state based on QC matrix
    inf = open(conf_dict['Step2_ExpMat']['qcmatfull'])
    for line in inf:
        if line.startswith('cellname'):
            continue
        ll = line.split()
        conf_dict['Mapping_stat']['umi_gene'] += int(ll[2])
        conf_dict['Mapping_stat']['cdsN'] += int(ll[3])
        conf_dict['Mapping_stat']['utr3N'] += int(ll[4])
        conf_dict['Mapping_stat']['utr5N'] += int(ll[5])
        conf_dict['Mapping_stat']['intronN'] += int(ll[6])
        conf_dict['Mapping_stat']['intergenicN'] += int(ll[7])
    inf.close()
    conf_dict['Mapping_stat']['totalreads'] = int(sp('wc -l %s'%(conf_dict['General']['barcode_reform']))[0].split()[0])    
    conf_dict['Mapping_stat']['q30reads'] = int(sp('wc -l %s'%(conf_dict['General']['bed']))[0].split()[0])

    
    ### create  QC dir and conduct QC
    wlog('generate reads QC measurement with own script, based on sample down reads',logfile)
    qcdir = conf_dict['General']['outputdirectory'] + 'QC/'
    createDIR(qcdir)
    os.chdir(qcdir)
    conf_dict['QCplots'] = {}
    conf_dict['QCplots']['map_summary'] = qcdir + conf_dict['General']['outname'] + '_map_summary.txt'
    mapsummary_doc = """genomic region(Category)\treads number
total reads\t%s
mappble reads\t%s 
total UMI count\t%s
CDS exon UMI count\t%s
3'UTR UMI count\t%s
5'UTR UMI count\t%s
intron UMI count\t%s
intergenic UMI count\t%s
"""%(str(conf_dict['Mapping_stat']['totalreads']),
     str(conf_dict['Mapping_stat']['q30reads']),
     str(conf_dict['Mapping_stat']['umi_gene']),
     str(conf_dict['Mapping_stat']['cdsN']),
     str(conf_dict['Mapping_stat']['utr3N']),
     str(conf_dict['Mapping_stat']['utr5N']),
     str(conf_dict['Mapping_stat']['intronN']),
     str(conf_dict['Mapping_stat']['intergenicN']))
    outf = open(conf_dict['QCplots']['map_summary'],'w') 
    outf.write(mapsummary_doc)
    outf.close()
    ## reads quality
    t= time.time()
    readsqc(conf_dict['General']['sampledownsam'],conf_dict['General']['outname'])
    wlog('generate bulk cell QC measurement with own script, based on sample down reads',logfile)

    cmd = "bedtools intersect -a %s -b %s -c > %s"%(conf_dict['General']['outputdirectory'] + 'annotation/'+ conf_dict['General']['outname'] + '_gene_anno_binexon.bed', conf_dict['General']['sampledownbed'],conf_dict['General']['outname']+'_sampledown_on_gbbin.bed' )
    rwlog(cmd,logfile)
    GBcover(conf_dict['General']['outname']+'_sampledown_on_gbbin.bed',conf_dict['General']['outname'])
    cmd = "%s %s %s"%('Rscript',conf_dict['rscript']+'readsbulkQC.r',conf_dict['General']['outname'])
    rwlog(cmd,logfile)
    
#       cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_qul'],conf_dict['General']['sam'],conf_dict['General']['outname'])
#       rwlog(cmd,logfile)
#       ## reads nucleotide composition
#       cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_nvc'],conf_dict['General']['sam'],conf_dict['General']['outname'])
#       rwlog(cmd,logfile)
#       ## reads GC content
#       cmd = "%s -i %s -o %s"%(conf_dict['Step3_QC']['read_gc'],conf_dict['General']['sam'],conf_dict['General']['outname'])
#       rwlog(cmd,logfile)
#       readsqctime = time.time() -t
#       wlog("time for readsqc: %s"%(readsqctime),logfile)
#       ## reads genebody coverage
#       t= time.time()
#
#       cmd = "%s -i %s -o %s -r %s"%(conf_dict['Step3_QC']['gb_cover'],conf_dict['General']['sam'],conf_dict['General']['outname'],conf_dict['General']['outputdirectory'] + 'annotation/'+conf_dict['General']['outname']+'_gene_anno_fullbed.bed')
#       rwlog(cmd,logfile)
#       bulkqctime = time.time() -t
#       wlog("time for bulkqc: %s"%(bulkqctime),logfile)
#       mvcmd1 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.qual.heatmap.pdf',qcdir + conf_dict['General']['outname'] + '_quality_heatmap.pdf')
#       mvcmd2 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.NVC_plot.pdf',qcdir + conf_dict['General']['outname'] + '_NVC.pdf')
#       mvcmd3 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.GC_plot.pdf',qcdir + conf_dict['General']['outname'] + '_GC.pdf')
#       mvcmd4 = "mv %s %s"%(qcdir + conf_dict['General']['outname'] + '.geneBodyCoverage.pdf',qcdir + conf_dict['General']['outname'] + '_GBcover.pdf')
#       rwlog(mvcmd1,logfile)
#       rwlog(mvcmd2,logfile)
#       rwlog(mvcmd3,logfile)
#       rwlog(mvcmd4,logfile)
#

    conf_dict['QCplots']['read_qul'] = qcdir + conf_dict['General']['outname'] + '_Figure1_quality_heatmap.pdf'
    conf_dict['QCplots']['read_nvc'] = qcdir + conf_dict['General']['outname'] + '_Figure2_NVC.pdf'
    conf_dict['QCplots']['read_gc'] = qcdir + conf_dict['General']['outname'] + '_Figure3_GC.pdf'
    conf_dict['QCplots']['gb_cover'] = qcdir + conf_dict['General']['outname'] + '_Figure4_GBcover.pdf'
    bulkqctime = time.time() -t
    wlog("time for bulkqc: %s"%(bulkqctime),logfile)
    
    ### individual cell QC
    wlog('generate individual cell QC measurement',logfile)
    t = time.time()
    conf_dict['QCplots']['duprate'] = qcdir + conf_dict['General']['outname'] + '_Figure5_duprate.pdf'
    conf_dict['QCplots']['covergn'] = qcdir + conf_dict['General']['outname'] + '_Figure8_coverGN.pdf'
    conf_dict['QCplots']['intronrate'] = qcdir + conf_dict['General']['outname'] + '_Figure9_intronrate.pdf'

    if conf_dict['General']['png_for_dot'] == 1:
        conf_dict['QCplots']['umicovergn'] = qcdir + conf_dict['General']['outname'] + '_Figure7_umi_coverGN.png'
        conf_dict['QCplots']['cumumiduprate'] = qcdir + conf_dict['General']['outname'] + '_Figure6_cumUMI_duprate.png'
    else:
        conf_dict['QCplots']['umicovergn'] = qcdir + conf_dict['General']['outname'] + '_Figure7_umi_coverGN.pdf'
        conf_dict['QCplots']['cumumiduprate'] = qcdir + conf_dict['General']['outname'] + '_Figure6_cumUMI_duprate.pdf'        
 
    conf_dict['Step2_ExpMat']['qcmatcc'] = qcdir + conf_dict['General']['outname'] + "_qcmat_clustercell.txt" 
    conf_dict['Step2_ExpMat']['expmatcc'] = qcdir + conf_dict['General']['outname'] + "_expmat_clustercell.txt" 
    conf_dict['results']['expmatcc'] = qcdir + conf_dict['General']['outname'] + "_expmat_clustercell.txt" 

    if int(conf_dict['Step3_QC']['select_cell_measure']) ==1:
        use_cutoff = conf_dict['Step3_QC']['covergncluster']
    elif int(conf_dict['Step3_QC']['select_cell_measure']) ==2:
        use_cutoff = conf_dict['Step3_QC']['topumicellnumber']
    else:
        ewlog('select_cell_measure value can only be 1 or 2, current value is %s'%(conf_dict['Step4_Analysis']['select_cell_measure']),logfile)

    cmd = "%s %s %s %s %s %s %s %s %s %s %s %s %s"%('Rscript',conf_dict['rscript']+'individual_qc.r',conf_dict['Step2_ExpMat']['qcmat'],conf_dict['Step2_ExpMat']['expmat'],conf_dict['General']['outname'],conf_dict['Step3_QC']['select_cell_measure'],use_cutoff,conf_dict['Step3_QC']['remove_low_dup_cell'],conf_dict['Step3_QC']['non_dup_cutoff'],conf_dict['Mapping_stat']['umi_gene'],conf_dict['Step2_ExpMat']['qcmatcc'],conf_dict['Step2_ExpMat']['expmatcc'],conf_dict['General']['png_for_dot'])
    rwlog(cmd,logfile)
    individualqctime = time.time() -t
    wlog("time for individualqc: %s"%(individualqctime),logfile)
    wlog("Step3 bulk and individual cell QC DONE",logfile)
    return conf_dict








