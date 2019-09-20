#!/bin/bash

echo $1
module load matlab
matlab  -nosplash -nodisplay -singleCompThread -r "load('$1'); cumulativeKill(blood, '$2'); exit" 
