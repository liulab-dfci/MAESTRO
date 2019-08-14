#!/usr/bin/env python
"""Description
Setup script for Dr.seq  -- QC and analysis pipeline for Drop-seq data
Copyright (c) 2015 Shengen Hu <tarelahu@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
import os
import sys
import subprocess
from distutils.core import setup, Extension



def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make 1>/dev/null 2>&1 ')
    sp('chmod 755 *')
    os.chdir(curdir)
    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1       
def main(): 
    has_R = check_R()
    if has_R == 0:
	    print("ERROR: Dr.seq requires R & Rscript under default PATH", file=sys.stderr)
	    sys.exit()
        
    has_bedtools = check_bedtools()
    print('Intalling Dr.seq, may take serval minutes')
    if has_bedtools == 0:
        compile_bedtools()
        setup(name="Drseqpipe",
              version="1.0.3",
              description="Drseq: Drop-seq QC and analysis pipeline",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://Tarela@bitbucket.org/tarela/drseq',
              package_dir={'Drseqpipe' : 'lib'},
              packages=['Drseqpipe'],
              package_data={'Drseqpipe': ['Config/Drseq_template.conf',
                                      'Rscript/analysis.r',
                                      'Rscript/individual_qc.r',
                                      'Rscript/readsbulkQC.r',
                                      'Rscript/post_analysis.r'
                                         ]},
              scripts=['bin/Drseq.py','refpackage/bedtools/bin/bedtools'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
          )
        print('bedtools is not detected under default PATH, bedtools is also installed')
        print('Installation of Dr.seq is DONE')
    
    else:
        setup(name="Drseqpipe",
              version="1.0",
              description="Drseq: Drop-seq QC and analysis pipeline",
              author='Shengen Hu',
              author_email='Tarelahu@gmail.com',
              url='https://Tarela@bitbucket.org/tarela/drseq',
              package_dir={'Drseqpipe' : 'lib'},
              packages=['Drseqpipe'],
              package_data={'Drseqpipe': ['Config/Drseq_template.conf',
                                      'Rscript/analysis.r',
                                      'Rscript/individual_qc.r',
                                      'Rscript/readsbulkQC.r',
                                      'Rscript/post_analysis.r'
                                         ]},
              scripts=['bin/Drseq.py'],
                        
              classifiers=[
            'Development Status :: version1.0 finish',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Artistic License',
            'Operating System :: POSIX',
            'Programming Language :: Python',
            'Topic :: pipeline',
            ],
              requires=[],
          )
        print('Installation of Dr.seq is DONE')


if __name__ == '__main__':
    main()

