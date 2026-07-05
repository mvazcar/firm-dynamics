% =========================================================================
% Per-period firm return (the "return function") for Hopenhayn & Rogerson
% (1993), in DOLLARS, as a function of the state (s, n_lag) and the
% employment choice n':
%
%   return_fn = p*s*(n')^theta  -  w*n'  -  p*cf  -  w*labor_adjustment(n',n_lag)
%             =  revenue      - wage bill - fixed cost - firing cost
%
% Conventions (all as in the paper): productivity s enters in LEVELS
% (s = exp(log s)); the fixed cost cf is in goods (dollar cost p*cf); the
% firing cost is in wages (dollar cost w*g). Vectorizes over s, n_lag, n'.
%
% The 2-D value-function iteration (vfn_2d) evaluates this over the whole
% (state x choice) grid; this file states the economics in one place.
% =========================================================================

function R = return_fn(s, n_lag, nprime, p, w, tau, params)
R = p.*s.*nprime.^params.theta ...
    - w.*nprime ...
    - p.*params.cf ...
    - w.*labor_adjustment(nprime, n_lag, tau) ;
