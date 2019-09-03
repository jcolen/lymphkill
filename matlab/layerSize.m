function [lsize] = layerSize(mask)

count = 0;
total = 0;
for i = 1:size(mask, 3)
    num = nnz(mask(:, :, i));
    if num > 0
        count = count + 1;
        total = total + num;
    end
end

lsize = total / count;

end