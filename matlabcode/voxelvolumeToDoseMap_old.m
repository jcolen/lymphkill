function [dmap] = voxelvolumeToDoseMap(infov, infod, volume, logic)

dmap = zeros(infod.Rows, infod.Columns, infod.NumberOfFrames);
dimd = size(dmap).';
dimv = infov.Dimensions.';
dd = [infod.PixelSpacing(1); infod.PixelSpacing(2); infod.SliceThickness];
if isempty(infod.SliceThickness)
    dd = [infod.PixelSpacing(1); infod.PixelSpacing(2); infov.PixelDimensions(3)];
end
dv = infov.PixelDimensions.';
cornerd = infod.ImagePositionPatient;
cornerv = infov.ImagePositionPatient;

for x = 1:dimd(1)
    for y = 1:dimd(2)
        for z = 1:dimd(3)
            pos = ([y; x; z] - 1).*dd + cornerd;
            %pos = ([x; y; z] - 1).*dd + cornerd;
            pos = (pos - cornerv)./dv + 1;
            if valid(pos, dimv)
                pos = ceil(pos);
                %dmap(x, y, z) = volume(pos(1), pos(2), pos(3));
                dmap(x, y, z) = volume(pos(2), pos(1), pos(3));
            end
        end
    end
end

if logic
    dmap = logical(dmap);
end

end

function [vld] = valid(pos, dim)
    %vld = all(pos > 0) && pos(1) < dim(1) && pos(2) < dim(2) && pos(3) < dim(3);
    vld = all(pos > 0) && pos(1) < dim(2) && pos(2) < dim(1) && pos(3) < dim(3);
end
