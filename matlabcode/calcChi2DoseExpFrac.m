function [chi2] = calcChi2DoseExpFrac(doseData, killAt05, regenRate, ignoreNegatives, printInfo)

chi2 = 0;
count = 0;

x1 = 5.0;
x2 = 0.5;
k1 = 0.992;

for i=1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    count = count + 1;
    total = sum(doseData(i).BinCounts);
    n = doseData(i).Fractions;
    b = n / (x2 - x1) * (log(1 - k1) / x1 - log(1 - killAt05) / x2);
    a = -log(1 - k1)/x1 - b * x1 / n;
    killed = 0;
    for j=1:numel(doseData(i).BinCounts)
        x = doseData(i).BinLowerBounds(j);
        killed = killed + doseData(i).BinCounts(j) * ...
            (1 - exp(-n * x * (a + b * x)));
    end
    percent = killed / total;
    
    regIndex = find(regenRate(:, 1) < doseData(i).PreTxLYA, 1, 'last');
    fpDays = doseData(i).Day - 25;
    if fpDays > 105
        fpDays = 105;
    end
    if fpDays > 0
        percent = percent - fpDays * regenRate(regIndex, 2);
    end
    
    diff = (percent - doseData(i).Measured) * doseData(i).PreTxLYA;
    if printInfo == 1
        fprintf('%s\t%f\n', doseData(i).Name, percent);
    end
    chi2 = chi2 + diff * diff;
end

chi2 = chi2 / (count - 1);

end