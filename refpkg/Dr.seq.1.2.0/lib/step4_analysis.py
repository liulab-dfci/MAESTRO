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
                                   createDIR)
# --------------------------
# main 
# --------------------------
def step4_analysis(conf_dict,logfile):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''   
    # start
    # create section for 
    t = time.time()
    wlog('Step4: analysis',logfile)
    wlog('dimentional reduction + clustering with own script, based on selected STAMP barcodes',logfile)
# Rscript analysis.r expmat outname coverGN highvarZ selectPCcutoff rdnumber maxKnum
    analysisdir = conf_dict['General']['outputdirectory'] + 'analysis/'
    createDIR(analysisdir)
    os.chdir(analysisdir)

    conf_dict['Step4_Analysis']['clusterresult'] = analysisdir + conf_dict['General']['outname']+'_cluster.txt'
    conf_dict['QCplots']['gapstat'] = analysisdir + conf_dict['General']['outname']+'_Figure10_GapStat.pdf'
    conf_dict['QCplots']['cluster'] = analysisdir + conf_dict['General']['outname']+'_Figure11_cluster.pdf'
    conf_dict['QCplots']['silhouette'] = analysisdir + conf_dict['General']['outname']+'_Figure12_silhouetteScore.pdf'
    conf_dict['QCplots']['umicolor'] = analysisdir + conf_dict['General']['outname']+'_Figure13_totalUMIcolored.pdf'
    conf_dict['QCplots']['itrcolor'] = analysisdir + conf_dict['General']['outname']+'_Figure14_intronRate_colored.pdf'
    conf_dict['results']['pctable'] = analysisdir + conf_dict['General']['outname']+'_pctable.txt'    
    conf_dict['results']['cortable'] = analysisdir + conf_dict['General']['outname']+'_correlation_table.txt' 
    conf_dict['results']['features'] = analysisdir + conf_dict['General']['outname']+'_features_clustercell.txt'
    
    cmd = "%s %s %s %s %s %s %s %s %s %s %s %s %s"%('Rscript',conf_dict['rscript']+'analysis.r',conf_dict['results']['expmatcc'],conf_dict['General']['outname'],conf_dict['Step4_Analysis']['highvarz'],conf_dict['Step4_Analysis']['selectpccumvar'],conf_dict['Step4_Analysis']['rdnumber'],conf_dict['Step4_Analysis']['maxknum'],conf_dict['Step4_Analysis']['pctable'],conf_dict['Step4_Analysis']['cortable'],conf_dict['Step4_Analysis']['clustering_method'],conf_dict['Step4_Analysis']['custom_k'],conf_dict['Step4_Analysis']['custom_d'])
    rwlog(cmd,logfile)
    cmd = '%s %s %s %s %s'%('Rscript',conf_dict['rscript']+'post_analysis.r',conf_dict['Step4_Analysis']['clusterresult'],conf_dict['Step2_ExpMat']['qcmatcc'],conf_dict['General']['outname'])
    rwlog(cmd,logfile)
    analysisqctime = time.time() - t
    wlog("time for analysis qc: %s"%(analysisqctime),logfile)
    wlog("Step4 analysis QC DONE",logfile)
    
    return conf_dict


