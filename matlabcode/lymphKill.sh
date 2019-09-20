#!/bin/bash

module load matlab
echo $1
if [ -z "$8" ]; then
    matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); [blood, blood_1frac] = lymphKill(masks, $3, $4, $5, '$1/$6'); save('$7', 'blood', 'blood_1frac'); exit;"
else
    matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); load('$1/${1}_dicom_info.mat'); fprintf('Rescaling masks\n'); masks = mask_fix(masks, infod, '$1/$6'); [blood, blood_1frac] = lymphKill(masks, $3, $4, $5, '$1/$6'); save('$7', 'blood', 'blood_1frac'); save('$1/${1}_masks_fixed.mat', 'masks'); exit;"
fi

./runDoseSummary.sh $7
