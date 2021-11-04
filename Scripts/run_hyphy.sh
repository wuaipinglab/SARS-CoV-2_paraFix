#!/usr/bin/env bash


HYPHY_DIR="Data/sampled_hyphy_results/"

for f in $(ls $HYPHY_DIR/*/*.nexus)
do
    if [[ ! -f ${f}.FUBAR.json ]]
    then
        echo $f
        cd $(dirname $f)
        hyphy fubar --alignment $(basename $f)
        cd -
    fi
done

HYPHY_DIR="Data/nextstrain_hyphy_results/"

for f in $(ls $HYPHY_DIR/*/*.nexus)
do
    if [[ ! -f ${f}.FUBAR.json ]]
    then
        echo $f
        cd $(dirname $f)
        hyphy fubar --alignment $(basename $f)
        cd -
    fi
done
