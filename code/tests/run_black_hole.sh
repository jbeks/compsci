#!/bin/bash
# run with "bash", not "sh"

# set simulation parameters

# time longer than longest simulation
# (used to create reference without black hole)
MAXT=1.6e10

# interpolation method
INTTYPE="hermite"

# used for example testing
DT=1080000
DIST=(1e9 5e9 10e9 15e9)
VEL=(20 240 460 680)
# used for actual testing
#DT=43200
#DIST=(1e9 2.5e9 5e9 7.5e9 10e9 12.5e9 15e9)
#VEL=(20 130 240 350 460 570 680)

# determine separator
if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

# set script and file names
INNAME="orbits.in"
INDIRNAME="input"
OUTDIRNAME="output"$SEP"out"
SCRIPTDIR="code"$SEP"tests"
SCRIPTNAME="black_hole.py"

# determine root directory
FILEDIR=$(dirname $(readlink -f $0))
ROOTDIR=$(dirname $(dirname $FILEDIR))

# determine paths to script, input and output directory
SCRIPT=$ROOTDIR$SEP$SCRIPTDIR$SEP$SCRIPTNAME
INFILE=$ROOTDIR$SEP$INDIRNAME$SEP$INNAME
OUTDIR=$ROOTDIR$SEP$OUTDIRNAME

# set scriptname for data shortening
SHORTDIR="code"$SEP"misc"
SHORTNAME="shorten_sim_data.py"
SHORT=$ROOTDIR$SEP$SHORTDIR$SEP$SHORTNAME

# set filename for baseline run
SIMDIR="code"$SEP"simulation"
SIMNAME="nbody.py"
SIM=$ROOTDIR$SEP$SIMDIR$SEP$SIMNAME
BASELINE="baseline_"$DT".out"

# run baseline
python -O $SIM < $INFILE $INTTYPE $MAXT $DT -o $DT > $OUTDIR$SEP$BASELINE".tmp"
python -O $SHORT < $OUTDIR$SEP$BASELINE".tmp" > $OUTDIR$SEP$BASELINE
rm -f $OUTDIR$SEP$BASELINE".tmp"

# for all distances and velocities, run black hole simulation
for dist in ${DIST[@]}
do
    for vel in ${VEL[@]}
    do
        # echo current distance and velocity
        echo "dist: "$dist" ; vel: "$vel
        OUTNAME="simdata_"$DT"_"$dist"_"$vel".out"
        # run black hole simulation
        python -O $SCRIPT < $INFILE $INTTYPE 0 $DT -o $DT \
            $dist $vel > $OUTDIR$SEP$OUTNAME".tmp"
        # shorten data
        python -O $SHORT < $OUTDIR$SEP$OUTNAME".tmp" > $OUTDIR$SEP$OUTNAME
        rm -f $OUTDIR$SEP$OUTNAME".tmp"
    done
done

