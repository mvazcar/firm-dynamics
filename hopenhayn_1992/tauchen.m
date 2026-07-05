% This function discretizes the AR(1) productivity process for the single-
% type Hopenhayn (1992) model using Tauchen (1986). The process is
%   y(t) = (1-rho)*mu + rho*y(t-1) + eps(t),   eps ~ N(0, sig^2).
% Returns the grid s (N x 1) and transition matrix Pi (N x N), Pi(j,k) =
% Pr(y' = s(k) | y = s(j)).
%
% Optional m sets the grid half-width in unconditional standard deviations
% (Tauchen's "q"); it defaults to 5. The grid spans mu +/- m*sqrt(sig^2/(1-rho^2)).

function [s, Pi] = tauchen(mu, rho, sig, N, m)

if nargin < 5 || isempty(m)
    m = 5 ;
end
s  = zeros(N,1) ;
Pi = zeros(N,N) ;

s(1) = mu - m*sqrt(sig^2/(1-rho^2)) ;
s(N) = mu + m*sqrt(sig^2/(1-rho^2)) ;

step = (s(N)-s(1))/(N-1) ;

for i = 2:(N-1)
   s(i) = s(i-1) + step ;
end

for j = 1:N
    for k = 1:N
        if k == 1
            Pi(j,k) = cdf_normal((s(1) - (1-rho)*mu - rho*s(j) + step/2) / sig) ;
        elseif k == N
            Pi(j,k) = 1 - cdf_normal((s(N) - (1-rho)*mu - rho*s(j) - step/2) / sig) ;
        else
            Pi(j,k) = cdf_normal((s(k) - (1-rho)*mu - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - (1-rho)*mu - rho*s(j) - step/2) / sig) ;
        end
    end
end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2)) ;
