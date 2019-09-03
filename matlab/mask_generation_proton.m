function [masks] = mask_generation_proton(contours, infov, infod, inputfile)
fid = fopen(inputfile);
data = textscan(fid, '%d %d %d %f', 'headerLines', 1);

otherInd = -1;

initialized = 0;
all = cell(1, 1);

for i = 1:numel(contours)
    if initialized
        all{1} = all{1} + voxelvolumeToDoseMap(infov, infod, contours(i).Segmentation, 1);
    else
        all{1} = voxelvolumeToDoseMap(infov, infod, contours(i).Segmentation, 1);
        initialized = 1;
    end
end

all = logical(all{1});

masks = struct;
for i = 1:numel(data{1})
    masks(i).Name = contours(data{1}(i)).ROIName;
    masks(i).Mask = voxelvolumeToDoseMap(infov, infod, contours(data{1}(i)).Segmentation, 1);
    masks(i).GV = data{2}(i);
    masks(i).Stationary = data{3}(i);
    masks(i).CardiacOutput = data{4}(i); 
    masks(i).LayerSize = layerSize(masks(i).Mask);
    if masks(i).CardiacOutput == -1
        otherInd = i;
    end

    all = logical(all - masks(i).Mask); 
end

%Now remove duplicated voxels in OtherOrgans
if otherInd ~= -1
    for i = 1:numel(masks)
        if i ~= otherInd
            masks(otherInd).Mask = masks(otherInd).Mask & ~masks(i).Mask;
        end
    end
    all = logical(all - masks(otherInd).Mask);
end

j = numel(masks) + 1;
masks(j).Name = 'Remaining';
masks(j).Mask = all;
masks(j).GV = 0;
masks(j).Stationary = 0;
masks(j).CardiacOutput = -2;
masks(j).LayerSize = layerSize(all);

end
