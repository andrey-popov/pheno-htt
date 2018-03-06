#!/bin/bash

# Convenience script to run MadGraph for multiple configuration files

set -e

for config in $*
do
    echo -e "\e[1;31mStart running ${config}\e[0m"
    mg5_aMC ${config}
    echo -e "\n\n"
done
