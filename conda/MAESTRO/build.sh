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
pip install --upgrade pip
pip install sinto

# install MAESTRO/R
$R CMD INSTALL
