function [errors] = calcDoseErrorsLinear(doseData, regenRate, killAt05, stdev, ignoreNegatives)
%Returns Nx3 array of format [value, upper bound, lower bound]

errors = zeros(numel(doseData), 3);

ki = zeros(1, 3);
ki(1) = killAt05;
ki(2) = killAt05 + stdev;
ki(3) = killAt05 - stdev;

for i=1:numel(doseData)
    if ignoreNegatives && (doseData(i).Measured < 0 || doseData(i).Measured >= 1)
        continue;
    end
    regIndex = find(regenRate(:, 1) < doseData(i).PreTxLYA, 1, 'last');
    fpDays = doseData(i).Day - 25;
    percent = zeros(1, 3);
    for j = 1:3
        percent(j) = calcKillLin(doseData(i), ki(j));

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

function [kill] = calcKillLin(row, killAt05)

% kill_0 = 0;
% kill_30 = 0.9;
% kill_50 = 0.992;

step = row.BinLowerBounds(2) - row.BinLowerBounds(1);

% rate_0 = (killAt05 - kill_0) / (0.5) * step;
% rate_05 = (kill_30 - killAt05) / (2.5) * step;
% rate_30 = (kill_50 - kill_30) / (2.0) * step;

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

rate = rate_0;

total = sum(row.BinCounts);
killed = 0;

for j = 1:numel(row.BinCounts)
    x = row.BinLowerBounds(j);
%     if x < 0.5
%         rate = rate + rate_0;
%     elseif x < 3.0
%         rate = rate + rate_05;
%     elseif x < 5.0
%         rate = rate + rate_30;
%     end
   if x < 0.5
        rate = rate + rate_0 * step;
    elseif x < 2.0
        rate = rate + rate_05 * step;
    %elseif x < 2.0
    %    rate = rate + rate_10 * step;
    elseif x < 3.0
        rate = rate + rate_20 * step;
    elseif x < 4.0
        rate = rate + rate_30 * step;
    elseif x < 5.0
        rate = rate + rate_40 * step;
    end
    killed = killed + row.BinCounts(j) * rate;
end

kill = killed / total;

end