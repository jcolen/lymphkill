function [errors] = calcDoseErrorsExp(doseData, regenRate, killAt05, stdev, ignoreNegatives)
%Returns Nx3 array of format [value, upper bound, lower bound]

a = zeros(1, 3);
b = zeros(1, 3);

[a(1), b(1)] = calcABExp(killAt05);
[a(2), b(2)] = calcABExp(killAt05 + stdev);
[a(3), b(3)] = calcABExp(killAt05 - stdev);

%Best fit Nakamura
a(:) = 0.3813;
b(:) = 0.0968;

errors = zeros(numel(doseData), 3);

for i=1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    regIndex = find(regenRate(:, 1) < doseData(i).PreTxLYA, 1, 'last');
    fpDays = doseData(i).Day - 25;
    percent = zeros(1, 3);
    for j = 1:3
        percent(j) = calcKillExp(doseData(i), a(j), b(j));

        if fpDays > 105
            fpDays = 105;
        end
        if fpDays > 0
            percent(j) = percent(j) - fpDays * regenRate(regIndex, 2);
        end
    end
    fprintf('%s\t%f\t%f\t%f\n', doseData(i).Name, percent(1), percent(2), percent(3))
    errors(i, :) = percent(:);
end

end

function [a, b] = calcABExp(killAt05)
    x1 = 5;
    x2 = 0.5;
    k1 = 0.992;
    %b = (0.2 * log(20.) + 2 * log(1 - killAt05)) / 4.5;
    %a = 0.2 * log(20.) - 5 * b;
    %b = 2/15 * log(10.) + 0.8 * log(1 - killAt05);
    %a = log(10.) / 3 - 3 * b;
    b = (log(1 - killAt05) - x2 / x1 * log(1 - k1)) / (x1 * x2 - x2 * x2);
    a = -1 / x1 * log(1 - k1) - b * x1;
end

function [kill] = calcKillExp(row, a, b)

total = sum(row.BinCounts);
killed = 0;

for j = 1:numel(row.BinCounts)
    x = row.BinLowerBounds(j);
    killed = killed + row.BinCounts(j) * ...
        (1 - exp(-a * x - b * x * x));
end

kill = killed / total;

end