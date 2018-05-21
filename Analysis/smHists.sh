#!/bin/bash

set -e

# Set mtt resolution to default value of 20% if not given
resolution=${1:-"0.2"}

sourceDir="../Showering/SM-tt"
mtt-hists ${sourceDir}/ttbar_[1-9]* -r $resolution &
mtt-hists ${sourceDir}/ttbar_FSR-up* -r $resolution &
mtt-hists ${sourceDir}/ttbar_FSR-down* -r $resolution &
mtt-hists ${sourceDir}/ttbar_mt-up* -r $resolution &
mtt-hists ${sourceDir}/ttbar_mt-down* -r $resolution &

wait

outputDir="hists"
mkdir -p $outputDir

hadd ${outputDir}/ttbar.root output/ttbar_[1-9]*
hadd ${outputDir}/ttbar_FSR-up.root output/ttbar_FSR-up*
hadd ${outputDir}/ttbar_FSR-down.root output/ttbar_FSR-down*
hadd ${outputDir}/ttbar_mt-up.root output/ttbar_mt-up*
hadd ${outputDir}/ttbar_mt-down.root output/ttbar_mt-down*

rm -r output
