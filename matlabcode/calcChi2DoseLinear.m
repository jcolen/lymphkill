function [chi2] = calcChi2DoseLinear(doseData, killAt05, regenRate, ignoreNegatives, printInfo)

chi2 = 0;

kill_0 = 0;
%kill_10 = 0.40;
kill_20 = 0.65;
kill_30 = 0.88;
kill_40 = 0.97;
kill_50 = 0.992;

rate_0 = (killAt05 - kill_0) / 0.5;
rate_05 = (kill_20 - killAt05) / 1.5;
%rate_10 = (kill_20 - kill_10) / 1.0;
rate_20 = (kill_30 - kill_20) / 1.0;
rate_30 = (kill_40 - kill_30) / 1.0;
rate_40 = (kill_50 - kill_40) / 1.0;

count = 0;

for i=1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    count = count + 1;
    total = sum(doseData(i).BinCounts);
    killed = 0;
    rate = rate_0;
    step = doseData(i).BinLowerBounds(2) - doseData(i).BinLowerBounds(1);
    for j=1:numel(doseData(i).BinCounts)
        x = doseData(i).BinLowerBounds(j);
        if x < 0.5
            rate = rate + rate_0 * step;
        elseif x < 2.0
            rate = rate + rate_05 * step;
        elseif x < 3.0
            rate = rate + rate_20 * step;
        elseif x < 4.0
            rate = rate + rate_30 * step;
        elseif x < 5.0
            rate = rate + rate_40 * step;
        end
        killed = killed + doseData(i).BinCounts(j) * rate;
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
        fprintf('%s\t%f\t%f\t%f\t%f\n', doseData(i).Name, doseData(i).Measured, percent, doseData(i).Measured * doseData(i).PreTxLYA, percent * doseData(i).PreTxLYA);
    end
    chi2 = chi2 + diff * diff;
end

chi2 = chi2 / (count - 1);

end