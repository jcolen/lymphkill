#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1/$2'); load('$1/$3'); maskToFile(masks, '$1/$1_mask_file.txt'); fid = fopen('$1/$1_organs_file.txt', 'w'); for i = 1:numel(contours) fprintf(fid, '%d %s\n', i, contours(i).ROIName); end; fclose(fid); exit;"
