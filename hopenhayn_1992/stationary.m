%##########################################################################
%--------------------------------------------------------------------------
% Stationary equilibrium of the single-type Hopenhayn (1992) firm-dynamics
% model, recoded with the output price p and wage w explicit.
%
% This is the lean, textbook version obtained from the HNS22 code by:
%   * keeping a SINGLE firm type (no low/high mixture, no omega);
%   * removing labor-force growth (no g, no 1+g anywhere);
%   * denominating the entry cost ce and operating cost cf in GOODS
%     (dollar costs p*ce and p*cf) rather than in labor;
%   * dropping the lifecycle / age / concentration block and the
%     transition (feed_data) machinery.
%
% Prices p, w and the costs ce, cf are INPUTS (w is solved from free entry
% by solve_wage; ce is hardcoded upstream). z = p/w is reported as a
% diagnostic.
%
% Input struct fields: alpha, beta, tol, maxiter, ns, svec, F, G,
%                      p, w, ce, cf, mstar.
% Output struct fields: p, w, z, ce, cf, entry_value_check, sstar_ind,
%                      sstar, mstar, mustar, nstar, N, avg_fsize,
%                      avg_stsize, startup_rate, exit_rate, fsize_dist.
%--------------------------------------------------------------------------
%##########################################################################

function f = stationary(params)

% Get parameters
cf    = params.cf ;         % operating cost in GOODS
ns    = params.ns ;
svec  = params.svec ;
F     = params.F ;          % single transition matrix (columns sum to one)
G     = params.G ;          % startup productivity distribution
p     = params.p ;          % output price (numeraire)
w     = params.w ;          % wage (solved from free entry)
ce    = params.ce ;         % entry cost in GOODS (hardcoded input)
mstar = params.mstar ;      % mass of entrants (scale normalization)

% -------------------------------------------------------------------------
% Value function (dollars) and exit threshold
% -------------------------------------------------------------------------
V = vfn(F, p, w, cf, params) ;
sstar_ind = find(V) ;
if isempty(sstar_ind)
    results.flag = 0 ;
    f = results ;
    return
else
    sstar_ind = sstar_ind(1) ;
end

% Diagnostic: expected dollar value of an entrant. Free entry requires this
% to equal the dollar entry cost p*ce; solve_wage picks w to enforce it.
entry_value_check = sum(V.*G) ;

% Degenerate-equilibrium guard: warn if no firm ever exits.
if sstar_ind == 1
    warning('stationary:noExit', ...
        ['No firm exits (exit threshold at grid index 1) at z = p/w = %.4g. ', ...
         'The equilibrium is likely degenerate; check the hardcoded ce.'], p/w) ;
end

% -------------------------------------------------------------------------
% Stationary firm-size distribution (no growth):  mu = (I - Fnew)^{-1} m*Gnew
% -------------------------------------------------------------------------
% Startup distribution and transition matrix, net of exit below threshold.
Gnew = G ;
Gnew(1:sstar_ind-1) = 0 ;

Fnew = F ;
Fnew(:,1:sstar_ind-1) = zeros(ns, sstar_ind-1) ;

I = eye(ns) ;
mustar = (I - Fnew) \ (mstar*Gnew) ;
mustar(1:sstar_ind-1) = 0 ;

% -------------------------------------------------------------------------
% Employment and aggregates. With cf and ce in GOODS, the fixed cost and the
% entry cost do NOT employ labor, so firm employment is the production labor
% n*(s) and total labor N is just production labor summed over firms.
% -------------------------------------------------------------------------
nstar = n(params, svec, p, w) ;
nstar(1:sstar_ind-1) = 0 ;

N = sum(nstar.*mustar) ;

M          = sum(mustar) ;                         % mass of incumbents
avg_fsize  = sum(nstar.*mustar)/M ;                % average firm size (employment)
avg_stsize = sum(nstar.*Gnew)/sum(Gnew) ;          % average entrant size

% In a no-growth stationary equilibrium entry equals exit, so the startup
% rate and the exit rate coincide. The exit flow is the mass of incumbents
% transiting into the below-threshold (exiting) states; by mass balance it
% equals the surviving-entrant inflow mstar*sum(Gnew). We compute both and
% they cross-check.
Fmu = Fnew*mustar ;                                % next-period inflow by state
startup_rate = mstar*sum(Gnew)/M ;
exit_rate    = sum(Fmu(1:sstar_ind-1))/M ;         % flow into exiting states / incumbents

fsize_dist = mustar/M ;                            % firm-size distribution

% -------------------------------------------------------------------------
% Pack results
% -------------------------------------------------------------------------
results.p  = p ;
results.w  = w ;
results.z  = p / w ;                 % relative (real) output price
results.ce = ce ;
results.cf = cf ;
results.entry_value_check = entry_value_check ;    % = p*ce at the solved w
results.sstar_ind = sstar_ind ;
results.sstar     = svec(sstar_ind) ;              % exit-threshold productivity
results.mstar   = mstar ;
results.mustar  = mustar ;
results.nstar   = nstar ;
results.N       = N ;
results.avg_fsize    = avg_fsize ;
results.avg_stsize   = avg_stsize ;
results.startup_rate = startup_rate ;
results.exit_rate    = exit_rate ;
results.fsize_dist   = fsize_dist ;

f = results ;
