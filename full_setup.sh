#!/bin/bash
# run with "bash", not "sh"

if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

ROOT=$(dirname $(readlink -f $0))
CODE=$ROOT$SEP"code"
SIM=$CODE$SEP"simulation"

bash $ROOT$SEP"convert_data.sh"
bash $SIM$SEP"build_cmod.sh"

