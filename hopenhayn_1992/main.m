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
%
% CALIBRATION: the parameters below are the Hopenhayn & Rogerson (1993)
% calibration (see ../hopenhayn_rogerson_1993), so this Hopenhayn (1992) model and the HR1993
% model in ../hopenhayn_rogerson_1993 share the SAME parameterization and grid. The only thing
% that differs between the two is the model itself: hopenhayn_1992 uses the exit-AFTER-
% shock timing V(s) = max(0, prof + beta*E[V]), whereas hopenhayn_rogerson_1993 uses HR's
% exit-BEFORE-shock timing. Structural params (beta, theta, rho, sigma) come
% straight from the paper; a, cf, ce and the entrant distribution were
% calibrated in ../hopenhayn_rogerson_1993 to HR's Table 1/2 targets at p = w = 1.
% =========================================================================
clc
clear
close all

set(0,'defaulttextInterpreter','latex')

% ---- Numerical / preference parameters (HR1993) ------------------------
beta    = 0.8    ; % 5-year discount factor
alpha   = 0.64   ; % production curvature theta (f = s*n^theta)
ns      = 100    ; % Number of points in the productivity grid
tol     = 1e-8   ; % Convergence tolerance
maxiter = 1000   ; % VFI cap (beta = 0.8 converges well within this)

% ---- Prices and scale --------------------------------------------------
p     = 1 ;   % Output price (numeraire). Held fixed; the wage w is solved for.
mstar = 1 ;   % Mass of entrants (scale normalization)

% ---- Hardcoded costs, in GOODS (HR1993 calibration, from ../hopenhayn_rogerson_1993) -------
CE = 14.9087 ;   % entry cost     (goods)
CF = 15.1537 ;   % operating cost (goods)

% ---- Productivity process: log s' = a + rho*log s + eps ----------------
rho       = 0.93 ;                   % persistence of log s
sigma_eps = (1-alpha)*sqrt(0.53) ;   % innovation std (= 0.2621), from HR targets
mu        = 0.8707 ;                 % mean of log s (a = (1-rho)*mu = 0.0609)
nmax      = 5000 ;                   % top grid state implies ~5000 employees
gridR     = 5.5 ;                    % log-grid range below the top
vfrac     = 0.74 ;                   % entrants: uniform over bottom 74% of grid

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

% Fixed productivity grid (uniform in log s), top state at n*(s_max) = nmax,
% and transition matrix oriented so that columns of F sum to one.
logSmax = (1-alpha)*log(nmax) - log(alpha) ;    % n*(s_max) = nmax at p = w = 1
svec = linspace(logSmax - gridR, logSmax, ns)' ;
F    = transpose(tauchen_fixed(svec, mu, rho, sigma_eps)) ;
params.svec = svec ;
params.F    = F ;

% Startup productivity distribution G: uniform over the bottom vfrac of grid
nlow = floor(vfrac*ns) ;
G = zeros(ns,1) ;
G(1:nlow) = 1/nlow ;
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
save hopenhayn_1992_benchmark
