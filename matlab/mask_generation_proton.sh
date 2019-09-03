#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); infov = dicom_read_header('$1/$3'); infod = dicominfo('$1/$4'); masks = mask_generation_proton(contours, infov, infod, '$1/${1}_mask_input.txt'); save('$1/$5', 'masks'); exit;"
