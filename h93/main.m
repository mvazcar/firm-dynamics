% =========================================================================
% HOPENHAYN & ROGERSON (1993, JPE) -- baseline economy with NO firing tax
% (tau = 0), at which HR1993 reduces to the frictionless single-type
% Hopenhayn (1992) model. This folder (h93) reuses that model engine
% (n.m, prof_fn.m, vfn.m, entry_residual.m, solve_wage.m, stationary.m,
% tauchen.m) and supplies the HR1993 calibration below. The paper is in
% h93.pdf.
%
% Parameter values follow the VFI Toolkit replication of HR1993
% (github.com/vfitoolkit/vfitoolkit-matlab-replication, HopenhaynRogerson1993):
%   beta  = 0.8            discount factor
%   alpha = 0.64           production curvature (HR call this theta)
%   cf    = 12             fixed cost of production (in goods)
%   ce    = 40             entry cost (in goods)
%   log-AR(1) productivity: log z' = a + rho*log z + eps,  eps ~ N(0,sigma_eps^2)
%     a = 0.078,  rho = 0.93,  sigma_eps = sqrt((1-rho)*0.53)
%   Tauchen discretization with n_z = 20 points and q = 4
%   entrants drawn uniformly over the bottom 65% of productivity states
%
% Productivity enters production in levels, p*z*n^alpha; here z = exp(svec)
% with svec the log-z grid, so prof_fn's p*exp(s)*n^alpha equals p*z*n^alpha.
%
% NORMALIZATION: HR1993 fix the wage (w = 1) and let the output price clear
% free entry. This model instead fixes p = 1 (numeraire) and solves free
% entry E[V(z0)] = p*ce for the wage w. Only z = p/w matters, so the real
% allocation is identical; the reported w equals 1/(HR equilibrium price).
% =========================================================================
clc
clear
close all

set(0,'defaulttextInterpreter','latex')

% ---- Numerical / preference parameters (HR1993) ------------------------
beta    = 0.8    ; % Discount factor
alpha   = 0.64   ; % Production curvature (theta in HR1993)
ns      = 20     ; % Number of productivity grid points (n_z)
q       = 4      ; % Tauchen grid half-width in std. dev.
tol     = 1e-8   ; % Convergence tolerance
maxiter = 1000   ; % VFI cap (beta = 0.8 converges well within this)

% ---- Prices and scale --------------------------------------------------
p     = 1 ;   % Output price (numeraire). Held fixed; the wage w is solved for.
mstar = 1 ;   % Mass of entrants (scale normalization)

% ---- Hardcoded costs, in GOODS (HR1993) --------------------------------
CE = 40 ;   % entry cost   (goods)
CF = 12 ;   % operating cost (goods)

% ---- Productivity process: log z' = a + rho*log z + eps ----------------
a_const    = 0.078 ;
rho        = 0.93 ;
sigma_logz = sqrt(0.53) ;                       % target unconditional std of log z
sigma_eps  = sqrt((1-rho)*sigma_logz^2) ;       % innovation std (as in VFI Toolkit)
mu         = a_const/(1-rho) ;                  % long-run mean of log z

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

% Discretize the log-AR(1) (Tauchen, n_z points, half-width q) and orient
% so that columns of F sum to one. svec is the LOG-z grid; prof_fn/n use
% exp(svec) = z, so production is p*z*n^alpha as in HR1993.
[svec, F] = tauchen(mu, rho, sigma_eps, ns, q) ;
F = transpose(F) ;
params.svec = svec ;
params.F    = F ;

% Entrant productivity distribution: uniform over the bottom 65% of states.
nlow = floor(0.65*ns) ;
G = zeros(ns,1) ;
G(1:nlow) = 1/nlow ;
params.G = G ;

% ---- Solve for the equilibrium wage, then the stationary equilibrium ---
w = solve_wage(params) ;
params.w = w ;

stateqbm = stationary(params) ;

% ---- Report ------------------------------------------------------------
fprintf('\n=== Hopenhayn-Rogerson (1993) calibration, tau = 0 (p fixed, solve for w) ===\n') ;
fprintf('  p (output price, numeraire) = %.15g\n', stateqbm.p) ;
fprintf('  w (wage, solved)            = %.15g\n', stateqbm.w) ;
fprintf('  z = p/w (relative price)    = %.15g\n', stateqbm.z) ;
fprintf('  ce (hardcoded, goods)       = %.15g\n', stateqbm.ce) ;
fprintf('  cf (hardcoded, goods)       = %.15g\n', stateqbm.cf) ;
fprintf('  E[V(z0)] at solved w        = %.15g  (free entry: should equal p*ce)\n', ...
        stateqbm.entry_value_check) ;
fprintf('  p*ce (dollar entry cost)    = %.15g\n', stateqbm.p*stateqbm.ce) ;
fprintf('  entry-condition residual    = %.3e\n', ...
        stateqbm.entry_value_check - stateqbm.p*stateqbm.ce) ;
fprintf('  --\n') ;
fprintf('  exit-threshold productivity (log z) = %.6g (grid index %d of %d)\n', ...
        stateqbm.sstar, stateqbm.sstar_ind, ns) ;
fprintf('  average firm size (employment)      = %.6g\n', stateqbm.avg_fsize) ;
fprintf('  average entrant size                = %.6g\n', stateqbm.avg_stsize) ;
fprintf('  startup rate                        = %.6g\n', stateqbm.startup_rate) ;
fprintf('  exit rate                           = %.6g\n', stateqbm.exit_rate) ;
fprintf('  total employment N (mstar = %g)      = %.6g\n\n', mstar, stateqbm.N) ;

%% Save output
save h93_benchmark
