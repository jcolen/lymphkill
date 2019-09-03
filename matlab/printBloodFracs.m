fprintf('Name\t<0.5\t0.5-0.6\t0.6-0.7\t0.7-0.8\t0.8-0.9\t0.9-1.0\t1.0-2.0\t2.0-3.0\t>3.0\n');

bounds = [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0];

for i = 1:numel(doseData)
    if doseData(i).Measured < 0.0 || doseData(i).Measured >= 1.0
        continue;
    end
    total = sum(doseData(i).BinCounts);
    fracs = zeros(size(bounds));
    blbs = round(doseData(i).BinLowerBounds, 3);
    for j=1:numel(bounds)-1
        ind0 = find(blbs == bounds(j));
        ind1 = find(blbs == bounds(j+1)) - 1;
        fracs(j) = sum(doseData(i).BinCounts(ind0:ind1)) / total;
    end
    ind0 = find(doseData(i).BinLowerBounds == bounds(end));
    fracs(end) = sum(doseData(i).BinCounts(ind0:end)) / total;
    fprintf('%s', doseData(i).Name);
    for j=1:numel(fracs)
        fprintf('\t%f', fracs(j));
    end
    %fprintf('\t%f', sum(fracs));
    fprintf('\n');
end