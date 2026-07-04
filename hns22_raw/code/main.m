% =========================================================================
% This program solves the benchmark model in Hopenhayn, Neira and
% Singhania (2022)
% =========================================================================
%% -------------------------------------------------------------------------
% Initialize parameters and other model primitives
% -------------------------------------------------------------------------
clc
clear
close all

set(0,'defaulttextInterpreter','latex') 

beta  = 1/1.04 ; % Discount rate
alpha = 0.64   ; % Labor share
g     = 0.01   ; % Trend in the labor force
ns    = 100    ; % Number of points in grid
tol   = 1e-8   ; % Specify tolerance

zstar = 1 ;   % Equilibrium output price 
mstar = 100 ; % Mass of entrants in steady state

% Create a structure that we want to transfer to functions
params.beta = beta ;
params.alpha = alpha ;
params.g = g ;
params.tol = tol ;
params.ns = ns ;
params.zstar = zstar ;
params.mstar = mstar ;

%% Labor Force data
datadir = '../data_summary_stats/' ;
fname = 'US_laborforce.csv' ;
lf_filename = strcat(datadir, fname) ;

%% 
% -------------------------------------------------------------------------
% Estimated Parameters
% -------------------------------------------------------------------------
% Targets from data
targets.entry_rate_78to83         = 0.1314*100; % 1978-83 Average Entry Rate
targets.entrant_size_78to83       = 5.96;       % 1978-83 Average Entrant Size
targets.conc_0yrold_78to83        = 0.0489*100; % 1978-83 Average Concentration of Entrants
targets.five_year_growth_78to83   = 0.7386*100; % 1978-83 Average five-year growth rate
targets.five_year_exit_78to83     = 0.5501*100; % 1978-83 Average five year exit rate
targets.tenKplusempshare_78to83   = 0.2600*100; % 1978-83 Average employment share of firms size greater than 10,000 employees
targets.onetofourfirmshare_78to83 = 0.5598*100; % 1978-83 Average firm share of firms size between 1 to 4 employees
targets.oneKplusempshare_78       = 0.4164*100; % 1978-83 Average employment share of firms size greater than 1000 employees
targets.agg_conc_78to83           = 0.5235*100; % 1978-83 Average employment share of firms size greater than 250 employees

% Set parameters
s0        = -4.344376541584754;   % Mean of startup distribution (\mu_g)
sigma0    = 1.331137767741511;    % Standard deviation of startup distribution (\sigma_G^2)
rho       = 0.984150757243253;    % Persistence  of AR1
sigma_eps = 0.245520815536363;    % Standard deviation of AR1 innovations
mu_l      = -2.431373086987380;   % Long-run AR1 mean of high type  
mu_h      = -1.436111629482697;   % Long-run AR1 mean of high type
cf_l      = 2.298950284374937;    % Fixed cost - intercept
cf_h      = 24.308026243791222;   % Fixed cost - intercept
omega_h   = 0.359430723369395;    % Share of high type firms

%% ------------------------------------------------------------------------

params.cf_h = cf_h ;
params.cf_l = cf_l ;

params.sigma0 = sigma0 ;
params.omega_h = omega_h;

[svec, F_l, F_h] = tauchen2(mu_l, mu_h, rho, sigma_eps, ns) ; 
F_l = transpose(F_l) ; % Transpose so that column sum to one
F_h = transpose(F_h) ; % Transpose so that column sum to one

% Add to params
params.svec = svec ;
params.F_l = F_l ;
params.F_h = F_h ;

% Startup productivity distribution
step = svec(2) - svec(1) ;
G = normcdf(svec + step/2, s0, sigma0) ;
G(2:end-1) = G(2:end-1) - G(1:end-2) ;
G(end) = 1 - sum(G(1:end-1)) ;
params.G = G ;

% Solve for stationary equilibrium in 1940, store results in stateqbm
stateqbm = stationary(params) ; 
ce = stateqbm.ce;

% Feed labor force data series to obtain moments, store results in bm
results = feed_data(params, stateqbm, lf_filename) ; 
bm = results.tab ;

% Store parameters
x_bench(1) = ce        ;    
x_bench(2) = s0        ;   
x_bench(3) = sigma0    ;  
x_bench(4) = rho       ; 
x_bench(5) = sigma_eps ;
x_bench(6) = mu_l      ;
x_bench(7) = mu_h      ;
x_bench(8) = cf_l      ;
x_bench(9) = cf_h      ;
x_bench(10) = omega_h  ;

Value = x_bench';

%% -------------------------------------------------------------------------
% The case with No Rise in labor force growth: feeds a constant 1978 lfg
% % -------------------------------------------------------------------------
fname = 'US_laborforce_norise.csv' ;
lf_filename_norise = strcat(datadir, fname) ;

params_norise = params ; 
stateqbm_norise = stationary(params_norise) ;

results_norise = feed_data(params, stateqbm, lf_filename_norise) ;
norise = results_norise.tab ; 

%% -------------------------------------------------------------------------
% The case with No Transition in labor force: feeds a constant 2014 lfg.
% -------------------------------------------------------------------------
fname = 'US_laborforce_notransition.csv' ;
lf_filename_notrans = strcat(datadir, fname) ;

params_notrans = params ; 
stateqbm_notrans = stationary(params_notrans) ;

results_notrans = feed_data(params, stateqbm, lf_filename_notrans) ;
notrans = results_notrans.tab ; 

%% Save output
save benchmark