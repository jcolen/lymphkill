function [ killAt05, stdev ] = calcDoseKillParams(doseData, regenRate, type, ignoreNegatives)


baseMask = zeros(1, numel(doseData));
baseMask(1:floor(0.8 * numel(doseData))) = 1;
baseMask = logical(baseMask);
N = 5;

killAt05 = 0;
vals = zeros(1, N);

for i = 1:N
    rows = baseMask(randperm(numel(baseMask)));
    if type == 1
        learnfun = @(x)calcChi2DoseExp(doseData(rows), x, regenRate, ignoreNegatives, 0);
    elseif type == 2
        learnfun = @(x)calcChi2DoseLinear(doseData(rows), x, regenRate, ignoreNegatives, 0);
    elseif type == 3
        learnfun = @(x)calcChi2DoseExpFrac(doseData(rows), x, regenRate, ignoreNegatives, 0);
    end
    killLearn = fminsearch(learnfun, 0.5);
    if type == 1
        fprintf('On training set:\t%f\n\n', calcChi2DoseExp(doseData(rows), killLearn, regenRate, ignoreNegatives, 1));
        fprintf('On test set:\t%f\n\n', calcChi2DoseExp(doseData(~rows), killLearn, regenRate, ignoreNegatives, 1));
    elseif type == 2
        fprintf('On training set:\t%f\n\n', calcChi2DoseLinear(doseData(rows), killLearn, regenRate, ignoreNegatives, 1));
        fprintf('On test set:\t%f\n\n', calcChi2DoseLinear(doseData(~rows), killLearn, regenRate, ignoreNegatives, 1));
    elseif type ==3
        fprintf('On training set:\t%f\n\n', calcChi2DoseExpFrac(doseData(rows), killLearn, regenRate, ignoreNegatives, 1));
        fprintf('On test set:\t%f\n\n', calcChi2DoseExpFrac(doseData(~rows), killLearn, regenRate, ignoreNegatives, 1));
    end
    killAt05 = killAt05 + killLearn;
    vals(i) = killLearn;
end

killAt05 = killAt05 / N;
vals
stdev = std(vals);
if type == 1
    fprintf('\nOn overall set:\t%f\n\n', calcChi2DoseExp(doseData, killAt05, regenRate, ignoreNegatives, 1));
elseif type == 2
    fprintf('\nOn overall set:\t%f\n\n', calcChi2DoseLinear(doseData, killAt05, regenRate, ignoreNegatives, 1));
elseif type == 3
    fprintf('\nOn overall set:\t%f\n\n', calcChi2DoseExpFrac(doseData, killAt05, regenRate, ignoreNegatives, 1));
end

end