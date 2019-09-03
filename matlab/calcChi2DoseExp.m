function [chi2] = calcChi2DoseExp(doseData, killAt05, regenRate, ignoreNegatives, printInfo)

chi2 = 0;

x1 = 5;
x2 = 0.5;
k1 = 0.992;
b = (log(1 - killAt05) - x2 / x1 * log(1 - k1)) / (x1 * x2 - x2 * x2);
a = -1 / x1 * log(1 - k1) - b * x1;

%Best fit Nakamura
%a = 0.3813;
%b = 0.0968;

count = 0;

for i=1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    count = count + 1;
    total = sum(doseData(i).BinCounts);
    killed = 0;
    for j=1:numel(doseData(i).BinCounts)
        x = doseData(i).BinLowerBounds(j);
        killed = killed + doseData(i).BinCounts(j) * ...
            (1 - exp(-a * x - b * x * x));
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
        %fprintf('%s\t%f\t%f\t%f\t%f\n', doseData(i).Name, doseData(i).Measured, percent, doseData(i).Measured * doseData(i).PreTxLYA, percent * doseData(i).PreTxLYA);
        fprintf('%s\t%f\n', doseData(i).Name, percent);
    end
    chi2 = chi2 + diff * diff;
end

chi2 = chi2 / (count - 1);

end