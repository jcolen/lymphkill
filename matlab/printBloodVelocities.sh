#!/bin/bash

module load matlab
echo $1
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); printBloodVelocities(masks, '$1/$6'); exit;"
