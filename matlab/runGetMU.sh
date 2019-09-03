#!/bin/bash

module load matlab
matlab -nojvm -nodisplay -nosplash -singleCompThread -r "RT_Plan_info = dicom_read_header('$1');getTotalMU(RT_Plan_info);save('$2', 'RT_Plan_info');exit"
