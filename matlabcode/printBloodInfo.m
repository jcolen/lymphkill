function [] = printBloodInfo(filename, killAt05)

load(filename);

calcBloodFracs(blood);

percent = calcPercentKillExp(blood, killAt05);

fprintf('\nOverall Percent Kill:\t%f\n', percent);

end
