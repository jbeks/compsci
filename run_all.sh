#!/bin/bash
# run with "bash", not "sh"

# currently uses parameters for quick example test
# change parameters in run_black_hole.sh
# and run_out_orbits.sh (both in ./code/tests/)
# (parameters used for report are commented out in those files)

# determine separator
if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

# determine dictionary paths
ROOT=$(dirname $(readlink -f $0))
CODE=$ROOT$SEP"code"
TESTS=$CODE$SEP"tests"
MISC=$CODE$SEP"misc"

# set script names
SETUPNAME="full_setup.sh"
BHNAME="run_black_hole.sh"
EVALNAME="run_out_orbits.sh"
HMNAME="heatmap.py"

# determine script paths
SETUP=$ROOT$SEP$SETUPNAME
BH=$TESTS$SEP$BHNAME
EVAL=$TESTS$SEP$EVALNAME
HM=$MISC$SEP$HMNAME

# run all scripts
#bash $SETUP
#bash $BH
#bash $EVAL
python $HM

