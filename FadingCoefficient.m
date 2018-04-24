function [ Coefficients ] = FadingCoefficient(mean, var, m, n)
    miu = log((mean^2)/sqrt(var + mean^2));
    sigma = sqrt(log(var/mean^2 + 1));
    Coefficients = lognrnd(miu, sigma, m ,n);
end

