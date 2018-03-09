#!/bin/bash

set -e

rm -r SM-tt*/logs

for syst in FSR-up FSR-down
do
    for name in $(ls "SM-tt_${syst}")
    do
        newName=$(echo $name | sed -E "s/ttbar_([0-9]+)/ttbar_${syst}_\1/")
        mv -i "SM-tt_${syst}/${name}" "SM-tt/${newName}"
    done
done

rm -r SM-tt_FSR-*
