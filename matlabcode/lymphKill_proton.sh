#!/bin/bash

module load matlab
echo $1
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); [blood, blood_1frac] = lymphKill_proton(masks, $3, $4, $5, '$1/$6'); save('$7', 'blood', 'blood_1frac'); exit;"

./runDoseSummary.sh $7
