import yaml
import sys
import json
import errno
import os

def check_files_exist(config):
	"""
	A function to check files/folders specified in th config.yaml file
	exist. 
	Also check if all the fastqs specified in the SAMPLE_JSON file
	exist. 
	"""
	with open(config, "r") as f:
		config = yaml.load(f, Loader = yaml.FullLoader)
		platform = config['platform']
		index_dir = config['genome']['mapindex']
		gtf = config['genome']['gtf']
		rsem = config['genome']['rsem']
		samples_json = config['SAMPLES_JSON']
		if not os.path.exists(index_dir):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), index_dir)
		if not os.path.isfile(gtf):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), gtf)
		if not os.path.isfile(samples_json):
			raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), samples_json)
		if platform == "Smartseq2":
			if not os.path.exists('rsem'):
				raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), rsem)

		# check all fastqs exist
		f2 = open(samples_json)
		samples = json.load(f2)
		f2.close()
		for sample in samples.keys():
			for read in samples[sample]:
				for fastq in samples[sample][read]:
					if not os.path.isfile(fastq):
						print('{fastq} does not exist'.format(fastq= fastq))
						raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fastq)
					else:
						pass


check_files_exist('config.yaml')



