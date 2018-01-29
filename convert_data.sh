#!/bin/bash
# run with "bash", not "sh"

if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

ROOT=$(dirname $(readlink -f $0))
DIR=$ROOT$SEP"ref"
FILES="$(find $DIR -type f -name "*_data.ref")"
SCRIPT="convert_data.py"

echo $DIR
for FILE in $FILES
do
    python $ROOT$SEP$SCRIPT $(basename $FILE "_data.ref")
done

