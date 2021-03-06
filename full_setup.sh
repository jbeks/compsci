#!/bin/bash
# run with "bash", not "sh"

# determine separator
if [[ $(uname -s) == Linux ]]
then
    SEP="/"
else
    SEP="\\"
fi

# determine paths to scripts
ROOT=$(dirname $(readlink -f $0))
CODE=$ROOT$SEP"code"
SIM=$CODE$SEP"simulation"

# execute setup bash scripts
bash $ROOT$SEP"convert_data.sh"
bash $SIM$SEP"build_cmod.sh"

