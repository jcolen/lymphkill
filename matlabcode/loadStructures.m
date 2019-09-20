function [] = loadStructures(rtssfile, imagedir, structfile, ctfile, dosefile, infofile)

rtssheader = dicominfo(rtssfile, 'UseVRHeuristic', false);
imageheaders = loadDicomImageInfo(imagedir, rtssheader.StudyInstanceUID);

names = fieldnames(rtssheader.StructureSetROISequence);
for i = 1:numel(names)
    field = getfield(rtssheader.StructureSetROISequence, names{i});
    fprintf('%s\n', field.ROIName);
end

fprintf('Loading in contours\n');
contours = readRTstructures(rtssheader, imageheaders);
save(structfile, 'contours');

infov = dicom_read_header(ctfile);
infod = dicominfo(dosefile);
save(infofile, 'infov', 'infod');

for i = 1:numel(contours)
    fprintf('%d %s\n', i, contours(i).ROIName);
end

end
