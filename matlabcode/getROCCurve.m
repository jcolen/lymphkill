function [X, Y, T, AUC] = getROCCurve(doseData, regenRate, classification, method, ignoreNegatives)
% doseData = patient dose information
% regenRate = replenishment information
% classification = 0 for non lymphopenic, 1 for lymphopenic
% method = 0 for arctan, 1 for exp, 2 for linear

[killAt05, stdev] = calcDoseKillParams(doseData, regenRate, method, ignoreNegatives);
if method == 1
    errors = calcDoseErrorsExp(doseData, regenRate, killAt05, stdev, ignoreNegatives);
elseif method == 2
    errors = calcDoseErrorsLinear(doseData, regenRate, killAt05, stdev, ignoreNegatives);
elseif method == 3
    errors = calcDoseErrorsLinear(doseData, regenRate, killAt05, stdev, ignoreNegatives);
end

preTx = [doseData(:).PreTxLYA]';
kill = [doseData(:).Measured]';
scores = preTx - preTx .* errors(:, 1);
pred = preTx - preTx .* kill; 
for i = 1:numel(doseData)
   fprintf('%s\t%f\t%f\t%d\n', doseData(i).Name, pred(i), scores(i), classification(i));
end
[X, Y, T, AUC] = perfcurve(classification, scores, logical(0));

end