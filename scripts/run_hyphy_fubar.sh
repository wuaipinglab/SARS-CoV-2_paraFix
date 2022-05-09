#!/usr/bin/env bash

HYPHY_DIR="output/nextstrain_hyphy_results/"

for f in $(ls $HYPHY_DIR/*/*.nexus)
do
    errlog=$(dirname $f)/errors.log
    if [[ ! -f ${f}.FUBAR.json ]] || [[ -f $errlog ]]
    then
        cd $(dirname $f)
        hyphy fubar --alignment $(basename $f)
        cd -
    fi
done
