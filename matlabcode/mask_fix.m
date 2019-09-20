function [masks] = mask_fix(masks, infod, filename)

fid = fopen(filename);
filenames = textscan(fid, '%s');
filenames = filenames{1};

infonew = dicominfo(filenames{1});

for i = 1:numel(masks)
    masks(i).Mask = voxelvolumeToDoseMap(infod, infonew, masks(i).Mask, 1);
    masks(i).LayerSize = layerSize(masks(i).Mask);
end

end