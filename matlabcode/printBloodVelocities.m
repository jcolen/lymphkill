function [] = printBloodVelocities(masks, filename)
%masks is a struct containing masks, names, and GV status for each organ
%filename is the name of the file containing the dose dicom files

fid = fopen(filename);
filenames = textscan(fid, '%s');
filenames = filenames{1};
fclose(fid);
num_beams = numel(filenames);

%For determining which voxels to crop out
used_voxels = cell(1, 1);
initialized = 0;

fprintf('Mask size: %d %d %d\n', size(masks(1).Mask, 1), size(masks(1).Mask, 2), size(masks(1).Mask, 3));

%Load in data
for i = 1:num_beams
    tmp = dicomread(filenames{i});
    if i == 1
        abc = squeeze(tmp);
        fprintf('Dose size: %d %d %d\n', size(abc, 1), size(abc, 2), size(abc, 3));
    end
    if initialized
        used_voxels{1} = used_voxels{1} + squeeze(tmp);
    else
        used_voxels{1} = squeeze(tmp);
        initialized = 1;
    end
end

used_voxels = logical(used_voxels{1});
nvoxels = nnz(used_voxels);
layer_size = floor(nvoxels / size(used_voxels, 3));
blood_voxels = 3 * nvoxels;

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
masked_voxels = zeros(size(used_voxels));
ind = -1;
for j = 1:numel(masks)
    if masks(j).CardiacOutput == -1
        ind = j;
    else
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
        fprintf('LINE: %20s\tTime Voxel: %f\tBlood Velocity: %f\tCO: %f\n', masks(j).Name, masks(j).TimeVoxel, 0.25 / masks(j).TimeVoxel, masks(j).CardiacOutput);
    end
    masked_voxels = masked_voxels + masks(j).Mask;
end

masked_voxels = logical(masked_voxels);
unmasked_voxels = logical(used_voxels - masked_voxels);
layers = 0;
for i = 1:size(unmasked_voxels, 3)
    if nnz(unmasked_voxels(:, :, i)) > 0
        layers = layers + 1;
    end
end
remainingLayerSize = nnz(unmasked_voxels) / layers;

if ind > 0
    masks(ind).CardiacOutput = remainingCardiacOutput * nnz(masks(ind).Mask) / ...
        (nnz(masks(ind).Mask) + nnz(unmasked_voxels));
    masks(ind).TimeVoxel = H2H * masks(ind).LayerSize / (masks(ind).CardiacOutput * blood_voxels);
    remainingCardiacOutput = remainingCardiacOutput - masks(ind).CardiacOutput;
    if masks(ind).TimeVoxel < mintvox
        %fprintf('\tOverwriting TimeVoxel for %s, was %f\n', ...
        %    masks(ind).Name, masks(ind).TimeVoxel);
        masks(ind).TimeVoxel = mintvox;
    end
    fprintf('LINE: %20s\tTime Voxel: %f\tBlood Velocity: %f\tCO: %f\n', masks(ind).Name, masks(ind).TimeVoxel, 0.25 / masks(j).TimeVoxel, masks(ind).CardiacOutput);
end

default_time_voxel = H2H * remainingLayerSize / (remainingCardiacOutput * blood_voxels);
if default_time_voxel < mintvox
    %fprintf('\tOverwriting blood velocity for Remaining, was %f\n', default_time_voxel);
    default_time_voxel = mintvox;
end
fprintf('LINE: %20s\tTime Voxel: %f\tBlood Velocity: %f\tCO: %f\n', 'Remaining', default_time_voxel, 0.25 / default_time_voxel, remainingCardiacOutput);

end
