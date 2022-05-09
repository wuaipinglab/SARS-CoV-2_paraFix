#!/usr/bin/env bash

HYPHY_DIR="output/nextstrain_hyphy_results/"

for f in $(ls $HYPHY_DIR/*/*.nexus)
do
    errlog=$(dirname $f)/errors.log
    if [[ ! -f ${f}.MEME.json ]] || [[ -f $errlog ]]
    then
        cd $(dirname $f)
        hyphy meme --alignment $(basename $f)
        cd -
    fi
done
