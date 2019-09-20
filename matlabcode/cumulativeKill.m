function [] = cumulativeKill(blood, plan_name)

blbs = [0:0.1:5];
mp = (blbs(2) - blbs(1)) / 2.0;
fracs = zeros(size(blbs));

total = numel(blood);

for j = numel(blbs)-1:-1:1
    fracs(j) = nnz(blood > blbs(j)) - sum(fracs(j+1:end));
end

fracs = fracs / total; 

%Best fit Nakamura
a = 0.3813;
b = 0.0968;

killfrac = zeros(size(blbs));
fprintf('Non-Normalized Kill Contributions\n');
for i = 1:numel(fracs)
    x = blbs(i) + mp;
    kp = 1 - exp(-a * x - b * x * x);
    killfrac(i) = fracs(i) * kp;
    if killfrac(i) > 0
        fprintf('Kill Contribution for [%.1f, %.1f] = %f\n', x - mp, x + mp, killfrac(i));
    end
end

f = figure('visible', 'off');
area(blbs, killfrac * 100. / sum(killfrac));
xlabel('Blood Dose (Gy)')
ylabel('Kill Contribution (%)')
xlim([0, 4.5]);
ti = sprintf('Average Kill Contribution for %s blood dose', plan_name);
title(ti);

pngname = sprintf('%s_cumkill.png', plan_name);
fprintf('PNG Filename:\t%s\n', pngname);
saveas(gcf, pngname);

end
