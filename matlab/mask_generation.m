function [masks] = mask_generation(contours, infov, infod, inputfile)
fid = fopen(inputfile);
data = textscan(fid, '%d %d %d %f', 'headerLines', 1);

otherInd = -1;

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
end

%Now remove duplicated voxels in OtherOrgans
if otherInd ~= -1
    for i = 1:numel(masks)
        if i ~= otherInd
            masks(otherInd).Mask = masks(otherInd).Mask & ~masks(i).Mask;
        end
    end
end

end
