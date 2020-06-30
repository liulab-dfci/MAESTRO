#!/bin/bash

# install MAESTRO/python
$PYTHON setup.py install

# install giggle
cd refpkg/giggle
make
cp bin/giggle $PREFIX/bin/
cd ../..

# install rabit
cd refpkg/Rabit
./configure --prefix=$PREFIX CFLAGS=-I${PREFIX}/include
make
make install
cd ../..

# install sinto, the pypi version is not useful, let's
# do it through git
git clone https://github.com/timoast/sinto
cd sinto
git checkout 32d8733be9ba79372001318174d3612dc73c28b0 # for version 0.7.1
$PYTHON setup.py install
cd ../

# there are two dependencies in R DESCRIPTION
# that can't be found in conda-forge or bioconda
# channel. They are grid and Gmisc
# 

# install MAESTRO/R
$R -e 'devtools::install(".", upgrade="never")'


