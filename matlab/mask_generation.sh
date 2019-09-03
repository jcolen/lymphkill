#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); load('$1/$3'); masks = mask_generation(contours, infov, infod, '$1/${1}_mask_input.txt'); save('$1/$4', 'masks'); exit;"
