#!/bin/bash

CURDIR=${PWD}
FILEDIR=$(dirname $(readlink -f $0))
cd $FILEDIR
python setup.py build_ext --inplace
cd $CURDIR

