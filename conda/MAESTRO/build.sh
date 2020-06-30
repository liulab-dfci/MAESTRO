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

# install sinto
#$PYTHON -m pip install --upgrade pip --no-deps
$PYTHON -m pip install sinto --no-deps

# there are two dependencies in R DESCRIPTION
# that can't be found in conda-forge or bioconda
# channel. They are grid and Gmisc
# 

# install MAESTRO/R
$R -e 'devtools::install(".", upgrade="never")'


