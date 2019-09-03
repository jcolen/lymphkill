function [chi2] = calcChi2Nakamura(alpha, beta)

measured = [1, 0.6, 0.35, 0.12, 0.03, 0.008];
chi2 = 0;

for i = 1:numel(measured)
    x = i - 1;
    surv = exp(-alpha * x - beta * x^2);
    chi2 = chi2 + (surv - measured(i))^2;
end

fprintf('%f %f - Chi2 = %f\n', alpha, beta, chi2);

end