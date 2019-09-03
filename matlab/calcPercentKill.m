function [percent] = calcPercentKill(blood)
%Calculate the total percent of blood cells killed given a blood dose
%vector

kill_0 = 0;
kill_05 = 0.1;
kill_20 = 0.5;
kill_30 = 0.9;

maxDose = max(blood(:));
step = 0.1;
prevDose = 0.0;
dose = step;

rate_0 = (kill_05 - kill_0) / (0.5 / step);
rate_05 = (kill_20 - kill_05) / (1.5 / step);
rate_20 = (kill_30 - kill_20) / (1.0 / step);
rate_30 = 0.9;

killed = 0;
rate = rate_0;

while dose < maxDose
    count = sum(blood < dose & blood > prevDose);
    %fprintf('%f %f %d %f\n', dose, rate, count, count * rate);
    
    if dose < 0.5
        rate = rate + rate_0;
    elseif dose < 2.0
        rate = rate + rate_05;
    elseif dose < 3.0
        rate = rate + rate_20;
    end
    
    killed = killed + count * rate;
    
    prevDose = dose;
    dose = dose + step;
end

percent = killed / numel(blood);
fprintf('%d %d %f\n', ceil(killed), numel(blood), percent);

end