function [kill] = calcPercentKillAtan(blood, killAt05)

y1 = tan(pi*(killAt05 - 0.5));
y2 = tan(2*pi/5);
a = 0.4 * (y2 - y1);
b = 1.2 * y1 - 0.2 * y2;

bins = [0:0.1:5];
counts = zeros(size(bins));
counts(end) = nnz(blood > bins(end));
for i = numel(bins)-1:-1:1
    counts(i) = nnz(blood > bins(i)) - sum(counts(i+1:end));
end

killed = 0;
for i = 1:numel(bins)
    killed = killed + counts(i) * ...
        (atan(a * bins(i) + b) / 3.14 + 0.5);
end

kill = killed / sum(counts);

end