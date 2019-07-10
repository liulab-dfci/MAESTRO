#!/usr/bin/python
import sys
from shutil import copyfile
from optparse import OptionParser

parser = OptionParser()

#metadata_file_name = 'data/TF_human_data_information.txt'
#curr_dir = 'data/TF_human'
#named_dir = 'data/named'
#name_map_file_name = 'cistrome_id_to_name_map.txt'


parser.add_option("-m",
                  "--metadata_file_name",
                  dest="metadata_file_name",
                  help="path to TF_human_data_information.txt file")

parser.add_option("-i",
                  dest="curr_dir",
                  help="Input data directory")

parser.add_option("-o",
                  dest="named_dir",
                  help="Output data directory")

parser.add_option("-n",
                  dest="name_map_file_name",
                  help="ID to name mapping file")

(options, args) = parser.parse_args()
if not options.metadata_file_name:
    parser.error('TF_human_data_information.txt file not given')
if not options.curr_dir:
    parser.error('Input data directory not given')
if not options.named_dir:
    parser.error('Output data directory not given')
if not options.name_map_file_name:
    parser.error('ID to name mapping file not given')

metadata = {}

for l in open(metadata_file_name, 'r'):
    A = l.rstrip().split('\t')
    metadata[A[7]] = { 'id' : A[0], \
                       'gsm' : A[1], \
                       'species' : A[2], \
                       'cell_line' : A[3], \
                       'cell_type' : A[4], \
                       'tissue' : A[5], \
                       'factor' : A[6], \
                       'file_name' : A[7] }

prefixes = {}

file_name_map = {}

for f in metadata.keys():
    prefix = '.'.join( [ metadata[f]['cell_line'], \
                         metadata[f]['cell_type'], \
                         metadata[f]['tissue'], \
                         metadata[f]['factor'] ] ) \
             .replace(' ', '_') \
             .replace('/', '-')

    if prefix not in prefixes:
        prefixes[prefix] = 0

    new_file_name = prefix + '.' + str(prefixes[prefix]) + '.bed'
    file_name_map[f] = new_file_name
    prefixes[prefix] = prefixes[prefix] + 1


name_map_file = open(name_map_file_name, 'w')

for f in file_name_map:
    src = curr_dir + '/' + f
    dst = named_dir + '/' + file_name_map[f]
    copyfile(src, dst)
    name_map_file.write(f + '\t' + file_name_map[f] + '\n')

name_map_file.close()

