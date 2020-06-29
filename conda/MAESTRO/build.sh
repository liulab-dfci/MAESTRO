#!/bin/bash

# install MAESTRO/python
$PYTHON setup.py install

# install Dr.seq

cd refpkg/Dr.seq.1.2.0
$PYTHON setup.py install
cd ../..

# install giggle
cd refpkg/giggle
make
cd ../..

# install rabit
cd refpkg/Rabit
./configure --prefix=$PREFIX
make
make install
cd ../..

# install sinto
pip install --upgrade pip
pip install sinto

# install MAESTRO/R
$R CMD INSTALL
