#!/bin/bash
# run with "bash", not "sh"

DT=216000
INTTYPE="rk4"
DIST=(1e9 1e10)
VEL=(20 80)

if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

INNAME="orbits.in"
INDIRNAME="input"
OUTDIRNAME="output"
SCRIPTDIR="code"$SEP"tests"
SCRIPTNAME="black_hole.py"

FILEDIR=$(readlink -f $0)
ROOTDIR=$(dirname $(dirname $(dirname $FILEDIR)))

SCRIPT=$ROOTDIR$SEP$SCRIPTDIR$SEP$SCRIPTNAME
INFILE=$ROOTDIR$SEP$INDIRNAME$SEP$INNAME
OUTDIR=$ROOTDIR$SEP$OUTDIRNAME

for dist in ${DIST[@]}
do
    for vel in ${VEL[@]}
    do
        OUTNAME="data_dt"$DT"_dist"$dist"_vel"$vel".out"
        python -O $SCRIPT < $INFILE $INTTYPE 0 $DT -o $DT \
            $dist $vel > $OUTDIR$SEP$OUTNAME
    done
done

