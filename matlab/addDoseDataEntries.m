function [doseData, doseData_1frac] = addDoseDataEntries(doseData, doseData_1frac, index, filename, name, measured, day, preLYA)
bins = [0:0.1:5];
bins_1frac = [0:0.02:1];
doseData(index).Name = name;
doseData_1frac(index).Name = name;

load(filename);

doseData(index).BinLowerBounds = bins;
doseData_1frac(index).BinLowerBounds = bins_1frac;

counts = zeros(size(bins));
counts_1frac = zeros(size(bins_1frac));

counts(end) = nnz(blood > bins(end));
counts_1frac(end) = nnz(blood_1frac > bins(end));

for j = numel(bins)-1:-1:1
    counts(j) = nnz(blood > bins(j)) - sum(counts(j+1:end));
end

for j = numel(bins_1frac)-1:-1:1
    counts_1frac(j) = nnz(blood_1frac > bins_1frac(j)) - sum(counts_1frac(j+1:end));
end

fprintf('%s\n', name);
calcBloodFracs(blood);
fprintf('After 1 fraction\n');
calcBloodFracs(blood_1frac);
fprintf('\n');

doseData(index).BinCounts = counts;
doseData(index).Measured = measured;
doseData(index).Day = day;
doseData(index).PreTxLYA = preLYA;

doseData_1frac(index).BinCounts = counts_1frac;
doseData_1frac(index).Measured = measured;
doseData_1frac(index).Day = day;
doseData_1frac(index).PreTxLYA = preLYA;

end
