function [chi2] = calcChi2DoseV2(doseData, killAt05, regenRate, printInfo)

chi2 = 0;

%y1 = tan(pi*(killAt05 - 0.5));
%y2 = tan(2*pi/5);
%a = 0.4 * (y2 - y1);
%b = 1.2 * y1 - 0.2 * y2;

syms a;
eqn = (cot(0.1 / a) - cot(1 / a)) / 3 - 2 * (cot ((1 - killAt05) / a) - cot(1 / a));
sola = vpasolve(eqn, 1);
aa = double(sola);
b = (cot(0.1 / aa) - cot(1 / aa)) / 3;
c = cot(1 / aa);
d = 1 - 3.1415926 * aa / 2;

for i=1:numel(doseData)
    if doseData(i).Measured < 0
        continue;
    end
    total = sum(doseData(i).BinCounts);
    killed = 0;
    for j=1:numel(doseData(i).BinCounts)
        killed = killed + doseData(i).BinCounts(j) * ...
            (aa * atan(b * doseData(i).BinLowerBounds(j) + c) + d);
            %(atan(a * doseData(i).BinLowerBounds(j) + b) / 3.14 + 0.5);
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
    
%     if regenRate > 0
%         fpDays = doseData(i).Day - 25;
%         if fpDays > 105
%             fpDays = 105;
%         end
%         if fpDays > 0
%             percent = percent - fpDays * regenRate;
%         end
%     end
    
    diff = percent - doseData(i).Measured;
    if printInfo == 1
        fprintf('%s\t%f\t%f\t%f\t%f\n', doseData(i).Name, doseData(i).Measured, percent, diff, diff*diff);
    end
    chi2 = chi2 + diff * diff;
end

chi2 = chi2 / (numel(doseData) - 1);

end