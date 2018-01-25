#!/bin/bash

CURDIR=${PWD}
FILEDIR="$(dirname "$0")"
cd $FILEDIR"/simulation"
python setup.py build_ext --inplace
cd $CURDIR
