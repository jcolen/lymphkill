function [doseData] = addDoseDataEntry(doseData, i, filename, name, measured, day, preLYA)
bins = [0:0.1:5];
doseData(i).Name = name;
load(filename);
doseData(i).BinLowerBounds = bins;
counts = zeros(size(bins));
if ~exist('blood_clinical')
    blood_clinical = blood;
end
counts(end) = nnz(blood_clinical > bins(end));
for j = numel(bins)-1:-1:1
    counts(j) = nnz(blood_clinical > bins(j)) - sum(counts(j+1:end));
end
fprintf('%s\n', name);
calcBloodFracs(blood_clinical);
doseData(i).BinCounts = counts;
doseData(i).Measured = measured;
doseData(i).Day = day;
doseData(i).PreTxLYA = preLYA;

end
