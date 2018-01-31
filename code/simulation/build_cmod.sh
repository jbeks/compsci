#!/bin/bash
# build c-module in directory of this file

CURDIR=${PWD}
FILEDIR=$(dirname $(readlink -f $0))
# move to this files directory
cd $FILEDIR
# build c module
python setup.py build_ext --inplace
# move back to original directory
cd $CURDIR

