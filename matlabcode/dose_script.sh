#!/bin/bash

module load matlab
echo '$1'
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "dose_script('$1', '$2', $3, $4, '$5'); exit"
