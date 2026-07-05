% Tauchen (1986) transition matrix on a PRE-SPECIFIED, uniformly (in logs)
% spaced grid. Used to replicate HR1993's construction, where the grid is
% pinned by the requirement that the top state implies ~5000 employees
% (rather than by the usual mean +/- q*sigma rule).
%
% AR(1):  logs' = (1-rho)*mu + rho*logs + eps,  eps ~ N(0, sigma^2).
% Inputs : logs (N x 1 uniform grid), mu (mean of log s), rho, sigma.
% Output : Pi (N x N), Pi(j,k) = Pr(logs' = logs(k) | logs = logs(j)).

function Pi = tauchen_fixed(logs, mu, rho, sigma)
N  = numel(logs) ;
logs = logs(:) ;
h  = logs(2) - logs(1) ;          % uniform step (in logs)
Pi = zeros(N,N) ;
cond_mean = (1-rho)*mu + rho*logs ;   % E[logs' | logs_j], N x 1

for j = 1:N
    m = cond_mean(j) ;
    for k = 1:N
        if k == 1
            Pi(j,k) = Phi((logs(1) - m + h/2)/sigma) ;
        elseif k == N
            Pi(j,k) = 1 - Phi((logs(N) - m - h/2)/sigma) ;
        else
            Pi(j,k) = Phi((logs(k) - m + h/2)/sigma) - ...
                      Phi((logs(k) - m - h/2)/sigma) ;
        end
    end
end
Pi = Pi ./ sum(Pi,2) ;            % guard: rows sum to one

function p = Phi(x)
    p = 0.5*erfc(-x/sqrt(2)) ;
