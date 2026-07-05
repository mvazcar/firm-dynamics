% This function discretizes the continuous AR(1) process in Hopenhayn, Neira, and Singhania (2022)
% by using the method proposed by Tauchen (1986). The AR(1) process takes the following form:
% y(t) = (1-rho)*mu_i + rho*y(t-1) + eps(t), where eps ~ N(0,sig^2) and i
% is firm type {low, high}. Joint grid is determined by mu1 (low mean in main code).

function [s, Pi1, Pi2] = tauchen2(mu1,mu2,rho,sig,N)

m       = 5 ;
s       = zeros(N,1);
Pi1      = zeros(N,N);
Pi2      = zeros(N,N);

s(1)    = mu1 - m*sqrt(sig^2/(1-rho^2));
s(N)    = mu1 + m*sqrt(sig^2/(1-rho^2));

step    = (s(N)-s(1))/(N-1);

for i=2:(N-1)
   s(i) = s(i-1) + step; 
end

for j = 1:N
    for k = 1:N
        if k == 1
            if ~isreal((s(1) - (1-rho)*mu1 - rho*s(j) + step/2) / sig)
            keyboard
            end
            Pi1(j,k) = cdf_normal((s(1) - (1-rho)*mu1 - rho*s(j) + step/2) / sig) ;
        elseif k == N
            Pi1(j,k) = 1 - cdf_normal((s(N) - (1-rho)*mu1 - rho*s(j) - step/2) / sig);
        else
            Pi1(j,k) = cdf_normal((s(k) - (1-rho)*mu1 - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - (1-rho)*mu1 - rho*s(j) - step/2) / sig);
        end
    end
end

for j = 1:N
    for k = 1:N
        if k == 1
            if ~isreal((s(1) - (1-rho)*mu2 - rho*s(j) + step/2) / sig)
            keyboard
            end
            Pi2(j,k) = cdf_normal((s(1) - (1-rho)*mu2 - rho*s(j) + step/2) / sig) ;
        elseif k == N
            Pi2(j,k) = 1 - cdf_normal((s(N) - (1-rho)*mu2 - rho*s(j) - step/2) / sig);
        else
            Pi2(j,k) = cdf_normal((s(k) - (1-rho)*mu2 - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - (1-rho)*mu2 - rho*s(j) - step/2) / sig);
        end
    end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));