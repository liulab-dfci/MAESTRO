#!/home1/wangchenfei/miniconda3/bin/python
# Time-stamp: <2016-04-29 Shengen Hu>
"""
 <Dr.seq: quality control and analysis for Drop-seq>
    Copyright (C) <2016>  <Shengen Hu>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string
import argparse
import subprocess

# -----------------------------------
# custom package
# -----------------------------------
import Drseqpipe

### tool function
from Drseqpipe.Utility      import          (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   readAnnotation,
                                   textformat,
                                   CMD
                                   )
### read and generate config file
from Drseqpipe.parse_config import (gen_conf,
                                   read_conf,
                                   make_conf
                                   )     


                                   
# -------------------
# main step
# -------------------
from Drseqpipe.step0_integrate_data   import step0_integrate_data
from Drseqpipe.step1_generate_matrix         import step1_generate_matrix
from Drseqpipe.step3_QC    import step3_QC
from Drseqpipe.step4_analysis import step4_analysis
from Drseqpipe.step5_summary import step5_summary

# ------------------------------------
# Misc functions
# ------------------------------------

    
#### read options

class ChiLinParser(argparse.ArgumentParser):
    '''
    Ex version for argparse(parameter parser) , add raise error function .
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args():
    '''
    Read parameter 
    '''
    description = "Dr.seq -- a quality control and analysis pipeline for droplet sequencing"
    parser = ChiLinParser(description = description)
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### generate config file
    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
                                             description = "Drseq config file generation. Usage: Drseq.py gen  -n config_name.conf")
    template_parser.add_argument("-n","--name", dest="config_name",required = True,help="name of your config file : config_name.conf")

    ### run config file
    pipe_parser = sub_parsers.add_parser("run", help = "run pipeline with a config file input",
                                         description = "Run Drseq pipeline with a config file input")
    pipe_parser.add_argument("-c","--config", required = True,
                             help = "specify the config file, -c config_name.conf" )
    pipe_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "The pipeline will over write output result if the output folder is already exist " )
    pipe_parser.add_argument("--clean",dest='Clean' , default=False, action='store_true',
                             help = "remove intermediate result generated during Dr.seq,default is NO" )
    ### simple mode
    simple_parser = sub_parsers.add_parser("simple", help = "run Drseq using simple mode",
                                         description = "(Run Drseq pipeline using simple mode/command line mode) Usage: Drseq.py -a barcode.fastq -b reads.fastq -n outname --maptool STAR -g mm10_refgenes.txt --mapindex /home/user/STAR_index")
    simple_parser.add_argument("-b","--barcode", dest = 'barcode',required = True,
                             help = "[required] barcode fastq file before any filtering step, only accept .fastq format" )
    simple_parser.add_argument("-r","--reads",dest='reads',required = True,
                             help = "[required] reads fastq file, accept raw fastq input or aligned sam format, format(fastq,sam) fixed by extension(.fastq or .sam)" )
    simple_parser.add_argument("-n","--name", dest="name",required = True,
                             help="[required] name of you config file and output dir, name only , no extension. The output files will be named like name.pdf, name.txt ... ")
    simple_parser.add_argument("--cellbarcodelength",dest='CBL' ,default='12', 
                             help = "specify the length of your cell barcode , default is 12, 12(cellbarcode) + 8(umi) = 20 (barcodefastq)" )
    simple_parser.add_argument("--umilength",dest='UMIL',  default='8', 
                             help = "specify the length of your UMI , default is 8, 12(cellbarcode) + 8(umi) = 20 (barcodefastq)" )
    simple_parser.add_argument("-g","--gene_annotation",dest='GA', required = False,
                             help = "[required absolute path if you didn't specific it in template config file] gene annotation file, the annotation file can be download from UCSC, full annotation text format(see documents for detail), or users can download gene annotation file in hg38 and mm10 version from our homepage" )
    simple_parser.add_argument("--maptool",dest='maptool' , choices = ("STAR", "bowtie2"),default="STAR", 
                             help = "choose mapping software for alignment, default is STAR, you can also choose bowtie2 as another option. mapping tool should corresponded to mapindex, STAR vs. STAR index, bowtie2 vs. bowtie2 index. Dr.seq will check your total memory if you choose STAR, becuse STAR consumes huge memory for mapping though fast. If you don't have 40G for total memory and choose STAR, Dr.seq will exit to protect your server/computer. you can also turn off the checkmem option to run STAR direclty" )
    simple_parser.add_argument("--checkmem",dest='checkmem' , choices = ("0", "1"), default="1",
                             help = "(default is 1 (on), only take effect when maptool = STAR) Dr.seq will check your total memory (if turned on, set 1, default) to make sure its greater 40G if you choose STAR as mapping tool. We don't suggest to run STAR on Mac. You can turn off (set 0) this function if you prefer STAR regardless of your memory (which may cause crash down of your computer) " )
    simple_parser.add_argument("--mapindex",dest='mapindex',required = False,
                             help = "[required if you didn't specific it in template config file] mapping index folder, there should be a mm10.star folder under mapindex folder if you use STAR to map to mm10 genome, mm10.bowtie2 folder if use bowtie2. for bowtie2, the index file should named like mm10.1.bt2(see documents for detail)" )
    simple_parser.add_argument("--thread",dest='P' ,default='8', 
                             help = "number of alignment threads to launch, ignored for sam input" )
    simple_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "specify the config file to create output folder , this cmd will rm existing result if set True ~!! " )
    simple_parser.add_argument("--clean",dest='Clean' , default=False, action='store_true',
                             help = "remove intermediate result generated during Dr.seq,default is No" )
    simple_parser.add_argument("--select_cell_measure",dest='select_cell_measure' , choices = ("1", "2"), default="1",
                             help = "Method to select STAMPs from cell_barcodes, choose from 1 or 2 (default is 1), 1: Cell_barcodes with more than 1000 genes covered are selected. 2: Top 1000 cell_harcodes with highest umi count will be selected" )
    simple_parser.add_argument("--remove_low_dup_cell",dest='remove_low_dup_cell' , choices = ("0", "1"), default="0",
                             help = "discard cell barcodes with low duplicate rate (< 0.1), set 1 to turn on this function, set 0 (default) to turn off. May not effective for samples with low sequencing depth" )

    args = parser.parse_args()
    ## generate config file template 
    if args.sub_command == "gen":
        gen_conf(args.config_name)
        sys.exit(0)

    ## run Dr.seq pipeline with config file input
    if args.sub_command == "run":
        if os.path.isfile(args.config):
            return args
        else:
            print('ERROR : -c input is not a config file\n')
            print(pipe_parser.print_help())
            sys.exit()
    
    ## run Dr.seq pipeline with a simple mode, input parameter in command line
    if args.sub_command == "simple":
        if args.name.endswith('.conf'):
            args.name = args.name[:-5]
        make_conf(args.barcode,args.reads,args.name,args.fover,args.CBL,args.UMIL,args.GA,args.P,args.mapindex,args.checkmem,args.maptool,args.select_cell_measure,args.remove_low_dup_cell)
        args.config = args.name + '.conf'
        return args

# ------------------------------------
# Main function
# ------------------------------------

def main():

    args = parse_args()
    conf_dict = read_conf(args.config)
    ### read raw path of output dir, the startdir will be used when the input file is not in absolute path
    conf_dict['General']['startdir'] = os.getcwd()+'/'

    ### check output name and dir from input parameter
    if conf_dict['General']['outname'] == "":
        print('your outname cannot be left blank,exit')
        sys.exit(1)
    if "." in conf_dict['General']['outname']:
        oldname = conf_dict['General']['outname']
        newname = oldname.replace(".","-")
        conf_dict['General']['outname'] = newname
        print('replace outname from %s to %s for latex summary'%(oldname,newname))
    if conf_dict['General']['outputdirectory'] == "":
        conf_dict['General']['outputdirectory'] = conf_dict['General']['outname']
        print('output directory is blank, use outname as directory name and set output directory in current folder')
    if "~" in conf_dict['General']['outname']:
        print('ERROR: ~ cannot appeared in outname, current outname is %s'%(conf_dict['General']['outname']))
        sys.exit(1)
    if "~" in conf_dict['General']['outputdirectory']:
        print('ERROR: require absolute path for outputdirectory')
        sys.exit(1)
    if not conf_dict['General']['outputdirectory'].endswith('/'):
        conf_dict['General']['outputdirectory'] += '/'
    if not conf_dict['General']['outputdirectory'].startswith('/'):
        conf_dict['General']['outputdirectory'] = conf_dict['General']['startdir'] + conf_dict['General']['outputdirectory']
    
    ### creat output dir
    if os.path.isfile(conf_dict['General']['outputdirectory'].rstrip("/")):
        print('ERROR: name of your output dir is exist as a file, cannot create a dir,Dr.seq exit')
        sys.exit(1)
    elif os.path.isdir(conf_dict['General']['outputdirectory']):
        if not args.fover:
            print('ERROR: name of your output dir is exist as a dir, Dr.seq exit because overwrite function is turned off, you can add -f parameter to turn on overwite function')
            sys.exit(1)
        else: 
            print('name of your output dir is exist as a dir, overwrite is turned on, write output result in existing dir')
    else:
        os.system("mkdir %s"%(conf_dict['General']['outputdirectory']))

    
    ### move to output dir
    os.chdir(conf_dict['General']['outputdirectory'])
    ## cp config file to output folder
    cmd = 'cp %s .'%(conf_dict['General']['startdir']+args.config)
    CMD(cmd)
    ### specify the main progress log file
    logfile = conf_dict['General']['outputdirectory']+'progress_log.txt'
    ### remove existing log file. 
    if os.path.isfile(logfile):
        CMD('rm %s'%logfile)
        
    ### Rscript location 
    #CONFIG_TEMPLATE = os.path.join(Drseq_pipe.__path__[0], "Config/Drseq_template.conf")
    conf_dict['rscript'] = os.path.join(Drseqpipe.__path__[0], "Rscript/")#'/mnt/Storage3/CR/Dropseq/drseq/Rscript/'
    conf_dict['clean'] = args.Clean
        
    ### main step for Dr.seq , see individual script for detail note.
    # preparing step, integrate parameter, prepare for following step
    t = time.time()
    step0_integrate_data(conf_dict,logfile)
    # main data processing step, including mapping, generate expression matrix and QC matrix which is used in next step
    step1_generate_matrix(conf_dict,logfile)
    step1time = time.time() -t
    wlog("running time for expression matrix generation: %s"%(step1time),logfile)
    # QC step, including bulk RNAseq QC(option), individual cell QC
    t = time.time()
    step3_QC(conf_dict,logfile)
    step3time = time.time()-t
    wlog("running time for QC: %s"%(step3time),logfile)
    # analysis step, including  select cell, filter high variance gene, pca + t-SNE dimentional reduction, k-means + Gap stat clustering
    t = time.time()
    step4_analysis(conf_dict,logfile)
    step4time = time.time() -t
    wlog("running time for clustering: %s"%(step4time),logfile)
    # summary step, integrate all QC figure and expression matrix, generate qC report with latex
    step5_summary(conf_dict,logfile)
    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt Dr.seq\n")
        sys.exit(1)

