function [] = calcBloodFracs(blood)

bounds = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0];

total = numel(blood);

for i = 1:numel(bounds) - 1
    fr = nnz(blood >= bounds(i) & blood < bounds(i+1)) / total;
    fprintf('Fraction %f < Dose < %f:\t%f\n', bounds(i), bounds(i+1), fr);
end

fr = nnz(blood >= bounds(end)) / total;
fprintf('Fraction Dose > %f:\t\t%f\n', bounds(end), fr);

end
