% =========================================================================
% Hopenhayn & Rogerson (1993) with firing costs -- FIRM BLOCK, 2-D model.
% Driver: solves the stationary equilibrium for tau = 0, 0.1, 0.2 and reports
% the firm-side rows of the paper's Table 3.
%
% State (s, n_lag) = (current productivity, lagged employment); choice n'.
% Bellman (paper p.7), with p fixed as the numeraire and w solved from entry:
%   W(s,n_lag) = max_{n'} { return_fn(s,n_lag,n')
%                           + beta*max( E_{s'|s} W(s',n'), -w*tau*n' ) }
% Model pieces: return_fn.m, labor_adjustment.m, vfn_2d.m, solve_wage_2d.m,
% stationary_2d.m, tauchen_fixed.m. At tau = 0 the state n_lag drops out and
% this collapses to the 1-D benchmark (../hopenhayn_rogerson_1993/calibrate_benchmark).
%
% Calibration is shared with the hopenhayn_1992/hopenhayn_rogerson_1993 benchmark (paper structural params;
% a, cf, ce, entrant distribution calibrated to HR's Table 1/2 targets).
% =========================================================================
clc ; clear ; close all ;

% ---- Calibration (from the benchmark) ----------------------------------
params.beta  = 0.8 ; params.theta = 0.64 ; params.rho = 0.93 ;
params.sigma = (1-params.theta)*sqrt(0.53) ;
params.mu    = 0.8707 ;  params.cf = 15.1537 ;  params.ce = 14.9087 ;
params.p     = 1 ;
params.vfrac = 0.74 ;

% ---- Grids and numerics ------------------------------------------------
params.nz    = 100 ;     % productivity grid points  (100/250 matches HR tightly;
params.na    = 250 ;     % employment grid points     50/150 runs faster for iterating)
params.nzR   = 5.5 ;     % log-s grid range below the top
params.nmax  = 5000 ;    % top productivity state implies n* = 5000 at p = w = 1
params.tol   = 1e-7 ; params.maxit = 2000 ;

% ---- Productivity grid + transition (same construction as benchmark) ---
logSmax = (1-params.theta)*log(params.nmax) - log(params.theta) ;
params.logs  = linspace(logSmax - params.nzR, logSmax, params.nz)' ;
params.s     = exp(params.logs) ;
params.Pi    = tauchen_fixed(params.logs, params.mu, params.rho, params.sigma) ;

% ---- Employment grid: 0 plus log-spaced up to a bit above nmax ---------
params.n     = [0, logspace(log10(1e-3), log10(1.2*params.nmax), params.na-1)]' ;
params.rev0  = params.s * (params.n'.^params.theta) ;       % nz x na, s(j)*n(l)^theta   (revenue/p)
params.fire0 = max(0, params.n - params.n') ;          % na x na, max(0, n_lag(i) - n'(l))

% ---- Entrant distribution over s: uniform bottom vfrac -----------------
nlow = floor(params.vfrac*params.nz) ;
params.v  = zeros(params.nz,1) ; params.v(1:nlow) = 1/nlow ;

% ---- Solve for each tau -------------------------------------------------
taus = [0, 0.1, 0.2] ;
R = cell(1,numel(taus)) ;
for it = 1:numel(taus)
    tau   = taus(it) ;
    w     = solve_wage_2d(tau, params) ;
    R{it} = stationary_2d(tau, w, params) ;
    R{it}.tau = tau ; R{it}.w = w ; R{it}.z = params.p/w ;
    fprintf('tau=%.2f  w=%.4f  z=%.4f  avg_size=%.2f  exit=%.3f  turnover=%.3f\n', ...
            tau, w, params.p/w, R{it}.avg_size, R{it}.exit_rate, R{it}.turnover) ;
end

% ---- Household closure (pins the scale; gives employment + welfare) -----
H = household(R, params, 0.6) ;   % A calibrated so N = 0.6 at tau = 0

% ---- Report Table 3 -----------------------------------------------------
b = R{1} ;
fprintf('\n===== HR1993 Table 3 (grid nz=%d, na=%d) =====\n', params.nz, params.na) ;
fprintf('  %-26s %8s %8s %8s | %s\n', '', 'tau=0', 'tau=.1', 'tau=.2', 'HR (0/.1/.2)') ;
row  = @(name, f, hr) fprintf('  %-26s %8.3f %8.3f %8.3f | %s\n', name, f(R{1}), f(R{2}), f(R{3}), hr) ;
rowH = @(name, x, hr)  fprintf('  %-26s %8.3f %8.3f %8.3f | %s\n', name, x(1), x(2), x(3), hr) ;
row('Price (z=p/w, rel.)',   @(r) r.z/b.z,            '1.000 / 1.026 / 1.048') ;
rowH('Consumption (output)',  H.cons_idx,             '100 / 97.5 / 95.4') ;
rowH('Average productivity',  100*H.prod/H.prod(1),   '100 / 99.2 / 97.9') ;
rowH('Total employment',      H.emp_idx,              '100 / 98.3 / 97.5') ;
rowH('Utility-adj. consumption', H.uac_idx,           '100 / 98.7 / 97.2') ;
row('Average firm size',     @(r) r.avg_size,         '61.2 / 61.8 / 65.1') ;
row('Layoff costs/wage bill', @(r) r.layoff_wb,       '0 / .026 / .044') ;
row('Job turnover rate',     @(r) r.turnover,         '.30 / .26 / .22') ;
row('Serial corr log n',     @(r) r.serial_n,         '.92 / .94 / .94') ;
row('Var of growth',         @(r) r.var_growth,       '.55 / .45 / .39') ;
row('Exit rate',             @(r) r.exit_rate,        '.39 /  -  /  - ') ;

save hr1993_2d_out
