% =========================================================================
% Verification for the Hopenhayn-Rogerson (1993) calibration (tau = 0) run
% in the single-type Hopenhayn model engine.
%
%   TEST 1: with p = 1 and the HR1993 parameters, solve for w and confirm
%           free entry clears (E[V(z0)] = p*ce), the exit threshold is
%           interior, and entry equals exit in steady state.
%
%   TEST 2 (homogeneity): rescale p by lambda != 1; the solved wage must
%           scale by lambda (z = p/w unchanged) and every physical quantity
%           (threshold, N, sizes, entry/exit rate) must be invariant.
%
% Run from the h93/ directory (all functions live there).
% =========================================================================
clc ; clear ; close all ;

tol_rel = 1e-6 ;

params = build_params() ;

% ---- TEST 1: base equilibrium, p = 1 -----------------------------------
params.p = 1 ;
params.w = solve_wage(params) ;
eq1      = stationary(params) ;

fprintf('\n============ TEST 1: HR1993 equilibrium (p=1) ============\n') ;
fprintf('  solved w = %.12g   z = p/w = %.12g\n', eq1.w, eq1.z) ;
fprintf('  exit threshold index = %d / %d   (interior expected)\n', eq1.sstar_ind, params.ns) ;
fprintf('  avg_fsize = %.6g   avg_stsize = %.6g   startup = %.6g   exit = %.6g\n\n', ...
        eq1.avg_fsize, eq1.avg_stsize, eq1.startup_rate, eq1.exit_rate) ;

pass = true ;
pass = chk('free entry E[V]=p*ce' , eq1.entry_value_check, eq1.p*eq1.ce   , tol_rel, pass) ;
pass = chk('threshold interior>1' , double(eq1.sstar_ind > 1)       , 1   , 0      , pass) ;
pass = chk('threshold interior<ns', double(eq1.sstar_ind < params.ns), 1   , 0      , pass) ;
pass = chk('entry = exit (SS)'    , eq1.startup_rate     , eq1.exit_rate  , tol_rel, pass) ;
pass = chk('N > 0'                , double(eq1.N > 0)               , 1   , 0      , pass) ;

% ---- TEST 2: homogeneity in (p,w) --------------------------------------
lambda    = 2 ;
params2   = params ;
params2.p = lambda ;
params2.w = solve_wage(params2) ;
eq2       = stationary(params2) ;

fprintf('============ TEST 2: homogeneity (p = %g) ============\n', lambda) ;
fprintf('  solved w = %.12g  (expected %g)   z = p/w = %.12g\n', eq2.w, lambda*eq1.w, eq2.z) ;

pass = chk('w scales by lambda'  , eq2.w              , lambda*eq1.w      , tol_rel, pass) ;
pass = chk('z invariant'         , eq2.z              , eq1.z             , tol_rel, pass) ;
pass = chk('free entry E[V]=p*ce', eq2.entry_value_check, eq2.p*eq2.ce    , tol_rel, pass) ;
pass = chk('threshold invariant' , eq2.sstar_ind      , eq1.sstar_ind     , 0      , pass) ;
pass = chk('N invariant'         , eq2.N              , eq1.N             , tol_rel, pass) ;
pass = chk('avg_fsize invariant' , eq2.avg_fsize      , eq1.avg_fsize     , tol_rel, pass) ;
pass = chk('startup invariant'   , eq2.startup_rate   , eq1.startup_rate  , tol_rel, pass) ;

fprintf('\n=====================================================\n') ;
if pass
    fprintf('ALL CHECKS PASSED.\n') ;
else
    fprintf('*** SOME CHECKS FAILED (see above). ***\n') ;
end
fprintf('=====================================================\n') ;


% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function pass = chk(name, got, ref, tol_rel, pass)
    if tol_rel == 0
        ok = (got == ref) ; err = abs(got - ref) ;
    else
        ok = abs(got - ref) <= tol_rel*max(1, abs(ref)) ;
        err = abs(got - ref)/max(1, abs(ref)) ;
    end
    status = 'PASS' ; if ~ok, status = 'FAIL' ; pass = false ; end
    fprintf('  [%s] %-22s got=%.10g  ref=%.10g  relerr=%.2e\n', status, name, got, ref, err) ;
end

function params = build_params()
    % HR1993 calibration (VFI Toolkit replication)
    params.beta    = 0.8 ;
    params.alpha   = 0.64 ;
    params.ns      = 20 ;
    params.tol     = 1e-8 ;
    params.maxiter = 1000 ;
    params.mstar   = 1 ;
    params.ce      = 40 ;   % goods
    params.cf      = 12 ;   % goods

    a_const    = 0.078 ;
    rho        = 0.93 ;
    sigma_logz = sqrt(0.53) ;
    sigma_eps  = sqrt((1-rho)*sigma_logz^2) ;
    mu         = a_const/(1-rho) ;
    q          = 4 ;

    [svec, F] = tauchen(mu, rho, sigma_eps, params.ns, q) ;
    params.svec = svec ;
    params.F    = transpose(F) ;

    nlow = floor(0.65*params.ns) ;
    G = zeros(params.ns,1) ;
    G(1:nlow) = 1/nlow ;
    params.G = G ;
end
