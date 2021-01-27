import argparse
import lisa
import time, os
import sys
def CommandLineParser():
    parser=argparse.ArgumentParser(description = "This is a description of input args")
    parser.add_argument("-S","--species", dest = "species", default = "GRCh38", help = "species of the data")
    parser.add_argument("-I","--input", dest = "input", default = "", help = "path of lisa data directory")
    return parser.parse_args()


parser = CommandLineParser()
if parser.species == "GRCh38":
    species = "hg38"
elif parser.species == "GRCm38":
    species = "mm10"
lisadir = parser.input


lisa_path = os.path.dirname(lisa.__file__)
if os.path.isdir(lisa_path):
        if os.path.exists(os.path.join(os.path.dirname(lisa.__file__), 'data', '%s_1000_2.0.h5' %(species))):
            print('lisa data already exists.')
        else:
            start_time = time.time()
            print('Start to install lisa data:',time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
            cmd = "lisa install %s multi %s" %(species, lisadir)
            os.system(cmd)
            end_time = time.time()
            print("End:", end_time-start_time)
else:
    print('lisa not found in current environment!')
