function [] = maskToFile(masks, filename)

fid = fopen(filename, 'w');

fprintf(fid, 'GV\tStatic\tName\n');

for j=1:numel(masks)
    fprintf(fid, '%d\t%d\t%s\n', masks(j).GV, masks(j).Stationary, masks(j).Name);
end

fclose(fid);

end
