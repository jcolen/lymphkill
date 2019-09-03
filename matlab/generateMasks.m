function [masks] = generateMasks(contours, infov, infod, indices, gv, thoracic)
masks = struct;
for i = 1:numel(indices)
    masks(i).Name = contours(indices(i)).ROIName;
    masks(i).Mask = voxelvolumeToDoseMap(infov, infod, contours(indices(i)).Segmentation, 1);
    if gv(i)
        masks(i).GV = 1;
    else
        masks(i).GV = 0;
    end
    masks(i).Stationary = 0;
    masks(i).LayerSize = layerSize(masks(i).Mask);
    masks(i).BloodVelocity = 5000 / 30 / (masks(i).LayerSize * infod.PixelSpacing(1) * infod.PixelSpacing(2) / 100.0) * 10;
    masks(i).TimeVoxel = infod.SliceThickness / masks(i).BloodVelocity;
    if isempty(infod.SliceThickness)
        masks(i).TimeVoxel = infov.PixelDimensions(3) / masks(i).BloodVelocity;
    end
end

masks(thoracic).Stationary = 1;

end
