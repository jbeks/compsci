#!/bin/bash
# run with "bash", not "sh"

DT=43200
INTTYPE="hermite"
DIST=(1e9 20e9 40e9 60e9 80e9 100e9)
VEL=(20 350 680 1010 1340 1670 2000)

MAXT=1.6e10

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

FILEDIR=$(dirname $(readlink -f $0))
ROOTDIR=$(dirname $(dirname $FILEDIR))

SCRIPT=$ROOTDIR$SEP$SCRIPTDIR$SEP$SCRIPTNAME
INFILE=$ROOTDIR$SEP$INDIRNAME$SEP$INNAME
OUTDIR=$ROOTDIR$SEP$OUTDIRNAME

SHORTDIR="code"$SEP"misc"
SHORTNAME="shorten_sim_data.py"
SHORT=$ROOTDIR$SEP$SHORTDIR$SEP$SHORTNAME

SIMDIR="code"$SEP"simulation"
SIMNAME="nbody.py"
SIM=$ROOTDIR$SEP$SIMDIR$SEP$SIMNAME
BASELINE="baseline_"$DT".out"

python -O $SIM < $INFILE $INTTYPE $MAXT $DT -o $DT > $OUTDIR$SEP$BASELINE".tmp"
python -O $SHORT < $OUTDIR$SEP$BASELINE".tmp" > $OUTDIR$SEP$BASELINE
rm -f $OUTDIR$SEP$BASELINE".tmp"

for dist in ${DIST[@]}
do
    for vel in ${VEL[@]}
    do
        echo "dist: "$dist" ; vel: "$vel
        OUTNAME="simdata_"$DT"_"$dist"_"$vel".out"
        python -O $SCRIPT < $INFILE $INTTYPE 0 $DT -o $DT \
            $dist $vel > $OUTDIR$SEP$OUTNAME".tmp"
        python -O $SHORT < $OUTDIR$SEP$OUTNAME".tmp" > $OUTDIR$SEP$OUTNAME
        rm -f $OUTDIR$SEP$OUTNAME".tmp"
    done
done

