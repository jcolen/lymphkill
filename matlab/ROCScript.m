[X1, Y1, T1, AUC1] = getROCCurve(doseData, RegenRate_var, postTxLYA < thresh, 1, ignoreNegatives);
[X2, Y2, T2, AUC2] = getROCCurve(doseData, RegenRate_var, postTxLYA < thresh, 2, ignoreNegatives);
[X3, Y3, T3, AUC3] = getROCCurve(doseData_1frac, RegenRate_var, postTxLYA < thresh, 3, ignoreNegatives);
plot(X1, Y1, X3, Y3, X2, Y2)
legend({'Exponential', 'Fractionated' 'Linear'}, 'Location', 'southeast');
xlabel('False Positive Rate');
ylabel('True Positive Rate');

ti = sprintf('ROC for Lymphopenia classification (LYA Threshold = %.2f)', thresh);
title(ti);