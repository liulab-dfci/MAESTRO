# -*- coding: utf-8 -*-
# @Author: Dongqing Sun
# @E-mail: Dongqingsun96@gmail.com
# @Date:   2020-02-28 22:23:48
# @Last Modified by:   Dongqing Sun
# @Last Modified time: 2020-03-05 17:13:48


import sys,os
import argparse as ap

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Merge microfluidic QC log files into singlecell.txt. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--log-dir", dest = "log_dir", type = str, default = "",  
        help = "Directory where QC log files are stored.")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO", 
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")

    return parser.parse_args()


def main():

    myparser = CommandLineParser()
    log_folder = myparser.log_dir
    directory = myparser.directory

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    out_file = os.path.join(directory, "singlecell.txt")

    output = {}
    for logfile in os.listdir(log_folder):
        if logfile.endswith('.mapping.log'):
            sample = logfile[:-12]
            # total,mapped,duplicate,uniq,mito,promoters,peaks = 0,0,0,0,0,0,0
            uniq,promoter = 0,0
            line_id = 0
            for line in open(log_folder + logfile).readlines():
                line = line.strip().split(' ')
                line_id += 1
                if line_id == 1:
                    uniq = line[0]
                if line_id == 2:
                    promoter = line[0]
                # if line_id == 1:
                #     total = line[0]
                # if line_id == 5:
                #     mapped = line[0]
                # if line_id == 4:
                #     duplicate = line[0]
                # if line_id == 14:
                #     uniq = line[0]
                # if line_id == 15:
                #     mito = line[0]
                # if line_id == 16:
                #     promoters = line[0]
                # if line_id == 17:
                #     peaks = line[0]
            output[sample] = [uniq,promoter]

    outf = open(out_file,'w')
    # outf.write('Cell\tUnique\tPromoter\n')
    for k in output:
        outf.write(k+'\t'+'\t'.join(output[k])+'\n')
    outf.close()

if __name__ == "__main__":
    main()
