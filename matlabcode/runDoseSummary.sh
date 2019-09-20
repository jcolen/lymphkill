#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "printBloodInfo('$1', 0.2044); exit;"
