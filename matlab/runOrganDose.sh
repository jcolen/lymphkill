#!/bin/bash

module load matlab

if [ -z "$4" ]; then
    matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1');check_organ_dose(masks, '$2', '$3');exit"
else
    matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2');load('$1/${1}_dicom_info.mat'); masks = mask_fix(masks, infod, '$1/$3'); fprintf('Masks fixed\n'); check_organ_dose(masks, '$1/$3', '$4');exit"
fi
