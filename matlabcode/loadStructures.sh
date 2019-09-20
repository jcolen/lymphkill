#!/bin/bash

module load matlab
#matlab -nojvm -nodisplay -nosplash -singleCompThread -r "rtssfile = '$1'; imagedir = '$2'; rtssheader = dicominfo(rtssfile); imageheaders = loadDicomImageInfo(imagedir, rtssheader.StudyInstanceUID); fprintf('Loading in contours\n'); contours = readRTstructures(rtssheader, imageheaders); save('$3', 'contours'); infov = dicom_read_header('$4'); infod = dicominfo('$5'); save('$6', 'infov', infod'); exit"

matlab -nojvm -nodisplay -nosplash -singleCompThread -r "loadStructures('$1', '$2', '$3', '$4', '$5', '$6'); exit"
