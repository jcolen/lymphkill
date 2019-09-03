function [errors] = calcDoseErrorsExpFrac(doseData, regenRate, killAt05, stdev, ignoreNegatives)

a = zeros(1, 3);
b = zeros(1, 3);
errors = zeros(numel(doseData), 3);

ab = zeros(1, 3);
count = 0;

for i = 1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    count = count + 1;
    regIndex = find(regenRate(:, 1) < doseData(i).PreTxLYA, 1, 'last');
    fpDays = doseData(i).Day - 25;
    n = doseData(i).Fractions;
    percent = zeros(1, 3);
    [a(1), b(1)] = calcABExpFrac(killAt05, n);
    [a(2), b(2)] = calcABExpFrac(killAt05 + stdev, n);
    [a(3), b(3)] = calcABExpFrac(killAt05 - stdev, n);
    ab = ab + a ./ b * n;
    
    for j = 1:3
        percent(j) = calcKillExpFrac(doseData(i), a(j), b(j), n);

        if fpDays > 105
            fpDays = 105;
        end
        if fpDays > 0
            percent(j) = percent(j) - fpDays * regenRate(regIndex, 2);
        end
    end
    fprintf('%s\t%f\t%f\t%f\n', doseData(i).Name, percent(1), percent(2), percent(3));
    errors(i, :) = percent(:);
end

fprintf('Average alpha/beta ratios:\t%f %f %f\n', ab(1) / count, ab(2) / count, ab(3) / count);

end

function [a, b] = calcABExpFrac(killAt05, n)
    x1 = 5;
    x2 = 0.5;
    k1 = 0.992;
    b = n / (x2 - x1) * (log(1 - k1) / x1 - log(1 - killAt05) / x2);
    a = -log(1 - k1)/x1 - b * x1 / n;
end

function [kill] = calcKillExpFrac(row, a, b, n)

total = sum(row.BinCounts);
killed = 0;

for j = 1:numel(row.BinCounts)
    x = row.BinLowerBounds(j);
    killed = killed + row.BinCounts(j) * ...
        (1 - exp(-n * x * (a + b * x)));
end

kill = killed / total;

end