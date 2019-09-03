function [blood, blood_1frac] = lymphKill_proton(masks, beam_on, gated, fracs, filename)
%masks is a struct containing masks, names, and GV status for each organ
%beam_on is the time on per beam
%gated is 1 if gated, 0 otherwise
%fracs is the number of fractions to use
%filename is the name of the file containing the dose dicom files

fid = fopen(filename);
filenames = textscan(fid, '%s');
filenames = filenames{1};
fclose(fid);

num_beams = numel(filenames);
info = cell(1, num_beams);
dose_masks = cell(numel(masks), num_beams); %Stores dose inside masked areas
frac_dose_masks = cell(numel(masks), num_beams);    %Stores total dose delivered to blood column in a fraction

fprintf('Mask size: %d %d %d\n', size(masks(1).Mask, 1), size(masks(1).Mask, 2), size(masks(1).Mask, 3));

%Load in data
for i = 1:num_beams
    info{i} = dicominfo(filenames{i});
    tmp = dicomread(filenames{i});
    if i == 1
        abc = squeeze(tmp);
        fprintf('Dose size: %d %d %d\n', size(abc, 1), size(abc, 2), size(abc, 3));
    end
    dose = info{i}.DoseGridScaling * squeeze(double(tmp)) / fracs / beam_on;
    for j = 1:numel(masks)
        dose_masks{j, i} = dose(masks(j).Mask);
    end
end

%Basic calculation for blood flow rate:
%Blood density = 5000 cm^3 / blood_voxels
%Blood_flow_rate = 5000 cm^3 / 30 s
%Voxel flow rate = Blood_flow_rate * Cardiac_output / (Layer_size *Blood_density)
%This simplifies to Vfl = CO * blood_voxels / (30 * Layer_size)
%In GVs, density is ~n times higher, Vfl = CO * blood_voxels / (n * 30 * Layer_size)
%So TimeVoxel = (n * 30 * Layer_size) / (CO * blood_voxels)

GV_density = 8; %GV density factor
H2H = 30.0;  %Heart to heart time in seconds
mintvox = 0.01; %Corresponds to max blood velocity of 25 cm/s

remainingCardiacOutput = 1.0;
OtherInd = -1;
RemainingInd = -2;
used_voxels = masks(1).Mask;
initialized = 0;
for j = 1:numel(masks)
    used_voxels = used_voxels + masks(j).Mask;
    if masks(j).CardiacOutput == -1
        OtherInd = j;
    elseif masks(j).CardiacOutput == -2
        RemainingInd = j;
    end
end

used_voxels = logical(used_voxels);
nvoxels = nnz(used_voxels);
blood_voxels = 3 * nvoxels;

for j = 1:numel(masks)
    if j ~= OtherInd && j ~= RemainingInd
        masks(j).TimeVoxel = H2H * masks(j).LayerSize / (masks(j).CardiacOutput * blood_voxels);
        if masks(j).GV  
            masks(j).TimeVoxel = masks(j).TimeVoxel * GV_density;
        else
            remainingCardiacOutput = remainingCardiacOutput - masks(j).CardiacOutput;
        end
        %Velocity capped at 20 cm/s
        if masks(j).TimeVoxel < mintvox
            %fprintf('\tOverwriting blood velocity for %s, was %f\n', ...
            %    masks(j).Name, masks(j).TimeVoxel);
            masks(j).TimeVoxel = mintvox;
        end
        fprintf('\t%s\t%f\tCO: %f\n', masks(j).Name, masks(j).TimeVoxel, masks(j).CardiacOutput);
    end
end

if OtherInd > 0
    masks(OtherInd).CardiacOutput = remainingCardiacOutput * nnz(masks(OtherInd).Mask) / ...
        (nnz(masks(OtherInd).Mask) + nnz(masks(RemainingInd).Mask));
    masks(OtherInd).TimeVoxel = H2H * masks(OtherInd).LayerSize / (masks(OtherInd).CardiacOutput * blood_voxels);
    remainingCardiacOutput = remainingCardiacOutput - masks(OtherInd).CardiacOutput;
    if masks(OtherInd).TimeVoxel < mintvox
        %fprintf('\tOverwriting TimeVoxel for %s, was %f\n', ...
        %    masks(OtherInd).Name, masks(OtherInd).TimeVoxel);
        masks(OtherInd).TimeVoxel = mintvox;
    end
    fprintf('\t%s\t%f\tCO: %f\n', masks(OtherInd).Name, masks(OtherInd).TimeVoxel, masks(OtherInd).CardiacOutput);
end

masks(RemainingInd).CardiacOutput = remainingCardiacOutput;
masks(RemainingInd).TimeVoxel = H2H * masks(RemainingInd).LayerSize / (masks(RemainingInd).CardiacOutput * blood_voxels);
if masks(RemainingInd).TimeVoxel < mintvox
    %fprintf('\tOverwriting blood velocity for Remaining, was %f\n', masks(RemainingInd).TimeVoxel);
    masks(RemainingInd).TimeVoxel = mintvox;
end
fprintf('\tRemaining\t%f\tCO: %f\n', masks(RemainingInd).TimeVoxel, remainingCardiacOutput);

for i = 1:num_beams
    for j = 1:numel(masks)
        if ~masks(j).Stationary
            dose_masks{j, i} = dose_masks{j, i} * masks(j).TimeVoxel;
        end
    end
end

%Fill in all doses to hit blood_voxels entries in the blood matrix
%This means either padding with zeros for non-GV organs
%or repeating to account for increased density for GV organs
pad_fac = 2 * GV_density + 2;
for i = 1:num_beams
    for j = 1:numel(masks)
        %GVs get ~8x blood density of non-GV organs
        if masks(j).GV
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:floor(blood_voxels / pad_fac)) = ...
                zeros(floor(blood_voxels / pad_fac) - mvoxels, 1, 'double');
            dose_masks{j, i} = repmat(dose_masks{j, i}, GV_density, 1);
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:blood_voxels) = ...
                zeros(blood_voxels - mvoxels, 1, 'double');
        else
            mvoxels = numel(dose_masks{j, i});
            dose_masks{j, i}(mvoxels+1:blood_voxels) = ...
                zeros(blood_voxels - mvoxels, 1, 'double');
        end
    end
end

if gated
    beam_on = beam_on * 3;
end

fprintf('Precalculating total beam doses\n');
for i = 1:num_beams
    fprintf('\tPrecalculating dose for Beam %d of %d\n', i, num_beams);
    for j = 1:numel(masks)
        fprintf('\t\t%s Dose\n', masks(j).Name);
        frac_dose_masks{j, i} = zeros(blood_voxels, 1);
        tvox = masks(j).TimeVoxel;
        lsize = floor(masks(j).LayerSize);
        if masks(j).Stationary
            frac_dose_masks{j, i} = dose_masks{j, i} * beam_on;
            continue;
        end
        t = 0;
        while t < beam_on
            if ~gated || (gated && mod(t, 4.0) < 1.333)
                frac_dose_masks{j, i} = frac_dose_masks{j, i} + dose_masks{j, i};
            end
            frac_dose_masks{j, i} = circshift(frac_dose_masks{j, i}, [-lsize, 1]);
            t = t + tvox;
        end
        %fprintf('\t\tMax:\t%f\n', max(frac_dose_masks{j, i}));
        %fprintf('\t\tMean:\t%f\n', mean(frac_dose_masks{j, i}(frac_dose_masks{j, i} > 0)));
    end
end

blood = zeros(blood_voxels, 1, 'double');

for day = 1:fracs
    fprintf('Beginning dose for Day %d\n', day);
    for i = 1:num_beams
        fprintf('\tApplying dose for Beam %d of %d\n', i, num_beams);
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
