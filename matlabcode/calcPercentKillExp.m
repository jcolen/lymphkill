function [percent] = calcPercentKillExp(blood, killAt05)

x1 = 5;
x2 = 0.5;
k1 = 0.992;
b = (log(1 - killAt05) - x2 / x1 * log(1 - k1)) / (x1 * x2 - x2 * x2);
a = -1 / x1 * log(1 - k1) - b * x1;

%Best fit Nakamura
%a = 0.3813;
%b = 0.0968;

bins = [0:0.1:5];
counts = zeros(size(bins));

for j = numel(bins)-1:-1:1
    counts(j) = nnz(blood > bins(j)) - sum(counts(j+1:end));
end

total = sum(counts);
killed = 0;
for j = 1:numel(counts)
    x = bins(j);
    killed = killed + counts(j) * (1 - exp(-a * x - b * x * x));
    %fprintf('Dose: %f\tCount: %d\tPercent: %f\n', x, counts(j), counts(j) / total);
end

percent = killed / total;


end
