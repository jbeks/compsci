#!/bin/bash
# run with "bash", not "sh"

# determine separator
if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

# determine paths to reference data
ROOT=$(dirname $(readlink -f $0))
DIR=$ROOT$SEP"ref"
FILES="$(find $DIR -type f -name "*_data.ref")"
SCRIPT="convert_data.py"

# create a data file for each reference
echo $DIR
for FILE in $FILES
do
    python $ROOT$SEP$SCRIPT $(basename $FILE "_data.ref")
done

