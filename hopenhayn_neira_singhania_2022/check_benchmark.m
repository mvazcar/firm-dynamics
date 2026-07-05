% =========================================================================
% Verification for the recoded (explicit p, w) HNS22 model.
%
%   TEST 1 (regression): with p = 1 and ce hardcoded to the benchmark value,
%           solve for w and confirm the recode reproduces the numbers the
%           ORIGINAL z-based code produced (w must come out to 1, z = 1).
%
%   TEST 2 (homogeneity): rescale the output price p by a factor lambda != 1
%           and re-solve. Because only z = p/w matters, the solved wage must
%           scale by the same lambda (z unchanged) and EVERY physical moment
%           -- including the aggregate labor share -- must be invariant. This
%           is the check that catches a stray wage factor hiding at w = 1.
%
% Run from anywhere; paths are resolved relative to this file.
% =========================================================================
clc ; clear ; close all ;

here    = fileparts(mfilename('fullpath')) ;
datadir = fullfile(here, 'data_summary_stats', filesep) ;   % data lives in this folder
addpath(here) ;                                             % model functions are here

tol_rel = 1e-6 ;   % relative tolerance for the regression checks

% ---- Reference values from the original z-based code (zstar = 1) ---------
REF = struct( ...
    'ce'            , 0.012815402758375272 , ...
    'N'             , 11.395399284055822   , ...
    'sstar_ind_l'   , 68                    , ...
    'sstar_ind_h'   , 74                    , ...
    'startup_rate'  , 0.10934545496514875  , ...
    'avg_fsize'     , 25.641718648693647   , ...
    'avg_stsize'    , 5.5469480389631292   , ...
    'avg_stsizeU'   , 0.015184065475327031 , ...
    'conc_age1'     , 0.05191035049679224  , ...
    'St6'           , 0.45677018650354134  , ...
    'bm1978_startup', 0.14824891472792379  , ...
    'bm1978_exit'   , 0.10994999167793053  , ...
    'bm1978_AFSt'   , 20.282153569051477   , ...
    'bm1978_lshare' , 0.76609604495465389  , ...
    'bm1978_c250'   , 0.52032266200013577  , ...
    'bm2014_startup', 0.087465224693520888 , ...
    'bm2014_exit'   , 0.091868694645193177 , ...
    'bm2014_AFSt'   , 28.6702456093434     , ...
    'bm2014_lshare' , 0.72427181220595593  , ...
    'bm2014_c250'   , 0.59698691480822419  , ...
    'norise_AFSt'   , 28.439951513245724   , ...
    'notrans_AFSt'  , 26.68184863264522    ) ;

% =========================================================================
% Build the model primitives (same as main.m)
% =========================================================================
params = build_params() ;
CE     = 0.012815402758375272 ;   % hardcoded benchmark entry cost

% ---- TEST 1: default experiment, p = 1, solve for w ---------------------
params.p  = 1 ;
params.ce = CE ;
params.w  = solve_wage(params) ;
eq1       = stationary(params) ;

r        = feed_data(params, eq1, fullfile(datadir,'US_laborforce.csv')) ;              bm1      = r.tab ;
r        = feed_data(params, eq1, fullfile(datadir,'US_laborforce_norise.csv')) ;       norise1  = r.tab ;
r        = feed_data(params, eq1, fullfile(datadir,'US_laborforce_notransition.csv')) ; notrans1 = r.tab ;

i78 = find(bm1.year == 1978) ;
i14 = find(bm1.year == 2014) ;
t1  = cumprod(1 + eq1.cond_growth_age(1:5)) ;

fprintf('\n============ TEST 1: benchmark regression (p=1) ============\n') ;
fprintf('  solved w = %.15g   z = p/w = %.15g\n', eq1.w, eq1.z) ;
fprintf('  free-entry residual E[V(s0)]-w*ce = %.3e\n\n', eq1.entry_value_check - eq1.w*eq1.ce) ;

pass = true ;
pass = chk('w == 1'            , eq1.w                 , 1                 , tol_rel, pass) ;
pass = chk('z == 1'            , eq1.z                 , 1                 , tol_rel, pass) ;
pass = chk('free entry E[V]=wce', eq1.entry_value_check, eq1.w*eq1.ce     , tol_rel, pass) ;
pass = chk('ce'               , eq1.ce                , REF.ce            , tol_rel, pass) ;
pass = chk('N'                , eq1.N                 , REF.N             , tol_rel, pass) ;
pass = chk('sstar_ind_l'      , eq1.sstar_ind_l       , REF.sstar_ind_l   , 0      , pass) ;
pass = chk('sstar_ind_h'      , eq1.sstar_ind_h       , REF.sstar_ind_h   , 0      , pass) ;
pass = chk('startup_rate'     , eq1.startup_rate      , REF.startup_rate  , tol_rel, pass) ;
pass = chk('avg_fsize'        , eq1.avg_fsize         , REF.avg_fsize     , tol_rel, pass) ;
pass = chk('avg_stsize'       , eq1.avg_stsize        , REF.avg_stsize    , tol_rel, pass) ;
pass = chk('avg_stsizeU'      , eq1.avg_stsizeU       , REF.avg_stsizeU   , tol_rel, pass) ;
pass = chk('conc_age(1)'      , eq1.conc_age(1)       , REF.conc_age1     , tol_rel, pass) ;
pass = chk('St(6)'            , eq1.St(6)             , REF.St6           , tol_rel, pass) ;
pass = chk('5yr growth'       , t1(end)-1             , 0.76685843795649999, tol_rel, pass) ;
pass = chk('bm1978 startup'   , bm1.startup_rate(i78) , REF.bm1978_startup, tol_rel, pass) ;
pass = chk('bm1978 exit'      , bm1.exit_rate(i78)    , REF.bm1978_exit   , tol_rel, pass) ;
pass = chk('bm1978 AFSt'      , bm1.AFSt(i78)         , REF.bm1978_AFSt   , tol_rel, pass) ;
pass = chk('bm1978 lshare'    , bm1.agglshare(i78)    , REF.bm1978_lshare , tol_rel, pass) ;
pass = chk('bm1978 c250'      , bm1.frac250plus_emp(i78), REF.bm1978_c250 , tol_rel, pass) ;
pass = chk('bm2014 startup'   , bm1.startup_rate(i14) , REF.bm2014_startup, tol_rel, pass) ;
pass = chk('bm2014 exit'      , bm1.exit_rate(i14)    , REF.bm2014_exit   , tol_rel, pass) ;
pass = chk('bm2014 AFSt'      , bm1.AFSt(i14)         , REF.bm2014_AFSt   , tol_rel, pass) ;
pass = chk('bm2014 lshare'    , bm1.agglshare(i14)    , REF.bm2014_lshare , tol_rel, pass) ;
pass = chk('bm2014 c250'      , bm1.frac250plus_emp(i14), REF.bm2014_c250 , tol_rel, pass) ;
pass = chk('norise AFSt(end)' , norise1.AFSt(end)     , REF.norise_AFSt   , tol_rel, pass) ;
pass = chk('notrans AFSt(end)', notrans1.AFSt(end)    , REF.notrans_AFSt  , tol_rel, pass) ;

% ---- TEST 2: homogeneity in (p,w) ---------------------------------------
lambda      = 2 ;
params2     = params ;
params2.p   = lambda ;                 % rescale the output price
params2.w   = solve_wage(params2) ;
eq2         = stationary(params2) ;
r           = feed_data(params2, eq2, fullfile(datadir,'US_laborforce.csv')) ; bm2 = r.tab ;

fprintf('\n============ TEST 2: homogeneity (p = %g) ============\n', lambda) ;
fprintf('  solved w = %.15g  (expected %g)   z = p/w = %.15g\n', eq2.w, lambda, eq2.z) ;

pass = chk('w scales by lambda', eq2.w              , lambda*eq1.w       , tol_rel, pass) ;
pass = chk('z invariant'       , eq2.z              , eq1.z              , tol_rel, pass) ;
pass = chk('free entry E[V]=wce', eq2.entry_value_check, eq2.w*eq2.ce     , tol_rel, pass) ;
pass = chk('N invariant'       , eq2.N              , eq1.N              , tol_rel, pass) ;
pass = chk('startup invariant' , eq2.startup_rate   , eq1.startup_rate   , tol_rel, pass) ;
pass = chk('avg_fsize invar.'  , eq2.avg_fsize      , eq1.avg_fsize      , tol_rel, pass) ;
pass = chk('conc_age1 invar.'  , eq2.conc_age(1)    , eq1.conc_age(1)    , tol_rel, pass) ;
pass = chk('lshare 1978 invar.', bm2.agglshare(i78) , bm1.agglshare(i78) , tol_rel, pass) ;
pass = chk('lshare 2014 invar.', bm2.agglshare(i14) , bm1.agglshare(i14) , tol_rel, pass) ;
pass = chk('AFSt 2014 invar.'  , bm2.AFSt(i14)      , bm1.AFSt(i14)      , tol_rel, pass) ;

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
        ok = (got == ref) ;
        err = abs(got - ref) ;
    else
        ok = abs(got - ref) <= tol_rel*max(1, abs(ref)) ;
        err = abs(got - ref)/max(1, abs(ref)) ;
    end
    status = 'PASS' ; if ~ok, status = 'FAIL' ; pass = false ; end
    fprintf('  [%s] %-22s got=%.12g  ref=%.12g  relerr=%.2e\n', ...
            status, name, got, ref, err) ;
end

function params = build_params()
    params.beta  = 1/1.04 ;
    params.alpha = 0.64 ;
    params.g     = 0.01 ;
    params.ns    = 100 ;
    params.tol   = 1e-8 ;
    params.maxiter = 100 ;   % HNS22 value-function-iteration cap (see README)
    params.mstar = 100 ;

    s0        = -4.344376541584754;
    sigma0    = 1.331137767741511;
    rho       = 0.984150757243253;
    sigma_eps = 0.245520815536363;
    mu_l      = -2.431373086987380;
    mu_h      = -1.436111629482697;
    params.cf_l    = 2.298950284374937;
    params.cf_h    = 24.308026243791222;
    params.omega_h = 0.359430723369395;
    params.sigma0  = sigma0 ;

    [svec, F_l, F_h] = tauchen2(mu_l, mu_h, rho, sigma_eps, params.ns) ;
    params.svec = svec ;
    params.F_l  = transpose(F_l) ;
    params.F_h  = transpose(F_h) ;

    step = svec(2) - svec(1) ;
    G = normcdf(svec + step/2, s0, sigma0) ;
    G(2:end-1) = G(2:end-1) - G(1:end-2) ;
    G(end) = 1 - sum(G(1:end-1)) ;
    params.G = G ;
end
