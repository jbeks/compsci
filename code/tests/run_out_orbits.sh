#!/bin/bash
# run with "bash", not "sh"

# set timestep for wich script is run
DT=43200

# determine separator
if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

# set script and file names
OUTDIRNAME="output"$SEP"out"
RESDIRNAME="output"$SEP"res"
SCRIPTDIR="code"$SEP"tests"
SCRIPTNAME="out_orbits.py"

# determine root directory
FILEDIR=$(dirname $(readlink -f $0))
ROOTDIR=$(dirname $(dirname $FILEDIR))

# determine paths to script, input directory and output directory
SCRIPT=$ROOTDIR$SEP$SCRIPTDIR$SEP$SCRIPTNAME
OUTDIR=$ROOTDIR$SEP$OUTDIRNAME
RESDIR=$ROOTDIR$SEP$RESDIRNAME

# set baseline file location
BASELINE=$OUTDIR$SEP"baseline_"$DT".out"

# run initialization
echo "Initializaton..."
python $SCRIPT $BASELINE -init

# run script for all .out files
OUTFILES=$OUTDIR$SEP"*.out"
for fname in ${OUTFILES[@]}
do
    # do not run baseline file
    if [[ $fname != $BASELINE ]]
    then
        echo $fname
        python $SCRIPT $fname -rd $RESDIR
    fi
done

