#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "load('$1');fprintf('%s\n', '$1'); for i=1:numel(contours) fprintf('%d %s\n', i, contours(i).ROIName); end;exit"
