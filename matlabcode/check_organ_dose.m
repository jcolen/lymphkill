function [] = check_organ_dose(masks, filename, outfile)

fid = fopen(filename);
filenames = textscan(fid, '%s');
filenames = filenames{1};
fclose(fid);

fout = fopen(outfile, 'w');

num_beams = numel(filenames);
tissue = zeros(size(masks(1)));

voxelsize = 0;

for i = 1:num_beams
    info = dicominfo(filenames{i});
    tmp = dicomread(filenames{i});
    tissue = tissue + double(squeeze(tmp) * info.DoseGridScaling);
    voxelsize = info.PixelSpacing(1) * info.PixelSpacing(2);
    if isempty(info.SliceThickness)
        fprintf('No slice thickness supplied, using value of 3 mm');
        voxelsize = voxelsize * 3;
    else
        voxelsize = voxelsize * info.SliceThickness;
    end
end

voxelsize = voxelsize / 1000; %Go from mm^3 to cm^3
dosechecks = [5, 10, 15, 20];

fprintf(fout, 'Organ,Max Dose (Gy),Mean Dose (Gy),Total Volume (cm^3),Integral Dose(Gy)');
for j = 1:numel(dosechecks)
    fprintf(fout, ',V%d (cm^3)', dosechecks(j));
end
fprintf(fout, '\n');

for i = 1:numel(masks)
    maxdose = max(tissue(masks(i).Mask));
    meandose = sum(tissue(masks(i).Mask)) / nnz(masks(i).Mask);
    volume = nnz(masks(i).Mask) * voxelsize;
    intdose = meandose * volume;
    fprintf(fout, '%20s,', masks(i).Name);
    fprintf(fout, '%f,', maxdose);
    fprintf(fout, '%f,', meandose);
    fprintf(fout, '%f,', volume);
    fprintf(fout, '%f', intdose);
    for j = 1:numel(dosechecks)
        dosevol = nnz(tissue(masks(i).Mask) >= dosechecks(j)) * voxelsize;
        fprintf(fout, ',%f', dosevol / volume);
    end
    fprintf(fout, '\n');
end

maxdose = max(tissue(:));
meandose = sum(tissue(:)) / nnz(tissue);
volume = nnz(tissue) * voxelsize;
intdose = meandose * volume;
fprintf(fout, 'ALL,');
fprintf(fout, '%f,', maxdose);
fprintf(fout, '%f,', meandose);
fprintf(fout, '%f,', volume);
fprintf(fout, '%f', intdose);
for j = 1:numel(dosechecks)
    dosevol = nnz(tissue(:) >= dosechecks(j)) * voxelsize;
    fprintf(fout, ',%f', dosevol / volume);
end
fprintf(fout, '\n');


end
