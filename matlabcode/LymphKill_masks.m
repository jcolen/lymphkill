function [blood, blood_1frac] = LymphKill_masks(masks, beam_on, gated, filename)
%masks is a struct containing masks, names, and GV status for each organ
%masks should also contain information such as the time voxel for each
%organ
%beam_on is 27.875 for RC static field and 95.4 for RS dynamic field

fid = fopen(filename);
filenames = textscan(fid, '%s');
filenames = filenames{1};
fclose(fid);

%time_voxel = 1/4;
time_voxel_gv = 1/40;
default_time_voxel = 1/4;
fracs = 5;

%total_bins = beam_on / time_voxel;
total_bins = beam_on / default_time_voxel;

num_beams = numel(filenames);
info = cell(1, num_beams);
dose = cell(1, num_beams);  %Stores dose outside of masked areas
dose_masks = cell(numel(masks), num_beams); %Stores dose inside masked areas
frac_dose = cell(1, num_beams); %Stores total dose delivered to blood column in a fraction
frac_dose_masks = cell(numel(masks), num_beams);    %Stores total dose delivered to blood column in a fraction

%For determining which voxels to crop out
used_voxels = cell(1, 1);
initialized = 0;

%For tracking layer sizes
total_mask = zeros(numel(masks), 1);
count_mask = zeros(numel(masks), 1);
layer_mask = zeros(numel(masks), 1);

%Load in data
for i = 1:num_beams
    info{i} = dicominfo(filenames{i});
    tmp = dicomread(filenames{i});
    if initialized
        used_voxels{1} = used_voxels{1} + squeeze(tmp);
    else
        used_voxels{1} = squeeze(tmp);
        initialized = 1;
    end
    
    dose{i} = info{i}.DoseGridScaling * squeeze(double(tmp)) / total_bins / fracs;
    for j = 1:numel(masks)
        dose_masks{j, i} = dose{i}(masks(j).Mask);
    end
    for j = 1:numel(masks)
        if masks(j).GV
            dose_masks{j, i} = dose_masks{j, i} * time_voxel_gv / default_time_voxel;
        end
        dose{i}(masks(j).Mask) = 0.0;
    end
    %dose{i} = info{i}.DoseGridScaling * squeeze(double(tmp)) / fracs / beam_on;
    %for j = 1:numel(masks)
    %    dose_masks{j, i} = dose{i}(masks(j).Mask) * masks(j).TimeVoxel;
    %end
    %for j = 1:numel(masks)
    %    dose{i}(masks(j).Mask) = 0.0;
    %end
end

used_voxels = logical(used_voxels{1});
for i = 1:num_beams
    dose{i} = dose{i}(used_voxels);
    %dose{i} = dose{i}(used_voxels) * default_time_voxel;
end
nvoxels = nnz(used_voxels);
blood_voxels = 3 * nvoxels;
%Fill in all doses to hit blood_voxels entries in the blood matrix
for i = 1:num_beams
    for j = 1:numel(masks)
        %GVs get increased densities
        if masks(j).GV
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:floor(blood_voxels / 18)) = ...
                zeros(floor(blood_voxels / 18) - mvoxels, 1, 'double');
            dose_masks{j, i} = repmat(dose_masks{j, i}, 8, 1);  %Increased density for GVs
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:blood_voxels) = ...
                zeros(blood_voxels - mvoxels, 1, 'double');
        else
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:blood_voxels) = ...
                zeros(blood_voxels - mvoxels, 1, 'double');
        end
    end
    
    dose{i} = reshape(dose{i}, [], 1);
    dose{i}(nvoxels+1:blood_voxels) = ...
        zeros(blood_voxels - nvoxels, 1, 'double');
end

%Determine average layer size in the z direction for each organ
for i = 1:size(masks(1).Mask, 3)
    for j = 1:numel(masks)
        if nnz(masks(j).Mask(:, :, i)) > 0
            count_mask(j) = count_mask(j) + 1;
            total_mask(j) = total_mask(j) + nnz(masks(j).Mask(:, :, i));
        end
    end
end

layer_mask = floor(total_mask ./ count_mask);
layer_size = floor(nvoxels / size(used_voxels, 3));
if gated
    beam_on = beam_on * 3;
end

fprintf('Precalculating total beam doses\n');
for i = 1:num_beams
    fprintf('\tPrecalculating dose for Beam %d of %d\n', i, num_beams);
    frac_dose{i} = zeros(blood_voxels, 1);
    fprintf('\t\tMain Dose\n');
    t = 0;
    while t < beam_on
        if ~gated || (gated && mod(t, 4.0) < 1.333)
            frac_dose{i} = frac_dose{i} + dose{i};
        end
        frac_dose{i} = circshift(frac_dose{i}, [-layer_size, 1]);
        %t = t + time_voxel;
        t = t + default_time_voxel;
    end
    fprintf('\t\tMax: %f\n', max(frac_dose{i}));
    fprintf('\t\tMean:\t%f\n', mean(frac_dose{i}(frac_dose{i} > 0)));
    for j = 1:numel(masks)
        fprintf('\t\t%s Dose\n', masks(j).Name);
        frac_dose_masks{j, i} = zeros(blood_voxels, 1);
        tvox = masks(j).TimeVoxel;
        if masks(j).Stationary
            frac_dose_masks{j, i} = dose_masks{j, i} * beam_on / tvox;
            continue;
        end
        t = 0;
        while t < beam_on
            if ~gated || (gated && mod(t, 4.0) < 1.333)
                frac_dose_masks{j, i} = frac_dose_masks{j, i} + dose_masks{j, i};
            end
            frac_dose_masks{j, i} = circshift(frac_dose_masks{j, i}, [-layer_mask(j), 1]);
            t = t + tvox;
        end
        fprintf('\t\tMax: %f\n', max(frac_dose_masks{j, i}));
        fprintf('\t\tMean:\t%f\n', mean(frac_dose_masks{j, i}(frac_dose_masks{j, i} > 0)));
    end
end

blood = zeros(blood_voxels, 1, 'double');

for day = 1:fracs
    fprintf('Beginning dose for Day %d\n', day);
    for i = 1:num_beams
        fprintf('\tApplying dose for Beam %d of %d\n', i, num_beams);
        blood = blood + frac_dose{i};
        blood = blood(randperm(blood_voxels));
        for j = 1:numel(masks)
            blood = blood + frac_dose_masks{j, i};
            blood = blood(randperm(blood_voxels));
        end
    end
    if day == 1
        blood_1frac = blood;
    end
end

end
