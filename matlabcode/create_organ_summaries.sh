#!/bin/bash

for d in *; do
    if [ -d $d ]
    then
        echo "${d}/${d}_clinical_filenames.txt"
        if [ -f ${d}/${d}_clinical_filenames.txt ] && [ -f ${d}/${d}_masks_fix.mat ]
        then
            ./runOrganDose.sh ${d}/${d}_masks_fix.mat ${d}/${d}_clinical_filenames.txt organ_summaries/${d}_summary.csv 
        elif [ -f ${d}/${d}_clinical_filenames.tx ] && [ -f ${d}/${d}_masks_fix.mat ]
        then
            ./runOrganDose.sh ${d}/${d}_masks_fix.mat ${d}/${d}_clinical_filenames.tx organ_summaries/${d}_summary.csv
        fi
    fi
done
