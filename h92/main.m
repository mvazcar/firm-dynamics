% =========================================================================
% Single-type Hopenhayn (1992) firm-dynamics model, with the output price p
% and the wage w explicit. Stripped-down version of the HNS22 recode:
%   * one firm type (no low/high mixture);
%   * no labor-force growth (g removed);
%   * entry cost ce and operating cost cf denominated in GOODS;
%   * stationary equilibrium only (no lifecycle profiles, no transition).
%
% DEFAULT EXPERIMENT:
%   * Fix the output price p = 1 (numeraire).
%   * Hardcode the entry cost ce and operating cost cf (in goods); see below.
%   * Solve the free-entry condition E[V(s0)] = p*ce for the wage w.
% =========================================================================
clc
clear
close all

set(0,'defaulttextInterpreter','latex')

% ---- Numerical / preference parameters ---------------------------------
beta    = 1/1.04 ; % Discount rate
alpha   = 0.64   ; % Curvature of the revenue function
ns      = 100    ; % Number of points in the productivity grid
tol     = 1e-8   ; % Convergence tolerance
maxiter = 100    ; % Value-function-iteration cap (see README note)

% ---- Prices and scale --------------------------------------------------
p     = 1 ;   % Output price (numeraire). Held fixed; the wage w is solved for.
mstar = 100 ; % Mass of entrants (scale normalization)

% ---- Hardcoded costs, in GOODS -----------------------------------------
% Entry and operating costs are denominated in units of the output good, so
% their dollar costs are p*CE and p*CF. (Numerical values carried over from
% the HNS22 low type; recalibrate as needed.)
CE = 0.012815402758375272 ;  % entry cost   (goods)
CF = 2.298950284374937    ;  % operating cost (goods)

% ---- Productivity process (single "low" type) --------------------------
s0        = -4.344376541584754;   % Mean of startup distribution
sigma0    =  1.331137767741511;   % Std. dev. of startup distribution
rho       =  0.984150757243253;   % AR(1) persistence
sigma_eps =  0.245520815536363;   % AR(1) innovation std. dev.
mu        = -2.431373086987380;   % AR(1) long-run mean (low type)

% ---- Assemble the parameter struct -------------------------------------
params.beta    = beta ;
params.alpha   = alpha ;
params.ns      = ns ;
params.tol     = tol ;
params.maxiter = maxiter ;
params.p       = p ;
params.ce      = CE ;
params.cf      = CF ;
params.mstar   = mstar ;

% Discretize the AR(1) and orient so that columns of F sum to one
[svec, F] = tauchen(mu, rho, sigma_eps, ns) ;
F = transpose(F) ;
params.svec = svec ;
params.F    = F ;

% Startup productivity distribution G (discretized normal, mean s0)
step = svec(2) - svec(1) ;
G = normcdf(svec + step/2, s0, sigma0) ;
G(2:end-1) = G(2:end-1) - G(1:end-2) ;
G(end) = 1 - sum(G(1:end-1)) ;
params.G = G ;

% ---- Solve for the equilibrium wage, then the stationary equilibrium ---
w = solve_wage(params) ;
params.w = w ;

stateqbm = stationary(params) ;

% ---- Report ------------------------------------------------------------
fprintf('\n=== Hopenhayn (1992): equilibrium (p fixed, solve for w) ===\n') ;
fprintf('  p (output price, numeraire) = %.15g\n', stateqbm.p) ;
fprintf('  w (wage, solved)            = %.15g\n', stateqbm.w) ;
fprintf('  z = p/w (relative price)    = %.15g\n', stateqbm.z) ;
fprintf('  ce (hardcoded, goods)       = %.15g\n', stateqbm.ce) ;
fprintf('  cf (hardcoded, goods)       = %.15g\n', stateqbm.cf) ;
fprintf('  E[V(s0)] at solved w        = %.15g  (free entry: should equal p*ce)\n', ...
        stateqbm.entry_value_check) ;
fprintf('  p*ce (dollar entry cost)    = %.15g\n', stateqbm.p*stateqbm.ce) ;
fprintf('  entry-condition residual    = %.3e\n', ...
        stateqbm.entry_value_check - stateqbm.p*stateqbm.ce) ;
fprintf('  --\n') ;
fprintf('  exit-threshold productivity = %.6g (grid index %d)\n', ...
        stateqbm.sstar, stateqbm.sstar_ind) ;
fprintf('  total employment N          = %.6g\n', stateqbm.N) ;
fprintf('  average firm size           = %.6g\n', stateqbm.avg_fsize) ;
fprintf('  average entrant size        = %.6g\n', stateqbm.avg_stsize) ;
fprintf('  startup rate                = %.6g\n', stateqbm.startup_rate) ;
fprintf('  exit rate                   = %.6g\n\n', stateqbm.exit_rate) ;

%% Save output
save h92_benchmark
