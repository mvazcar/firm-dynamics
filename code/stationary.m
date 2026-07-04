%##########################################################################
%--------------------------------------------------------------------------
% This function takes a structure of model parameters as input, solves for
% the stationary equilibrium in Hopenhayn, Neira and Singhania (2022) and
% returns the outcome.
%
% RECODE NOTE (explicit p and w):
%   The original code carried a single composite price z = p/w (called
%   zstar) and, holding z = 1 fixed, COMPUTED the entry cost ce from the
%   value functions. This version instead takes the output price p, the
%   wage w and the entry cost ce as INPUTS (ce is hardcoded upstream; w is
%   solved from free entry by solve_wage). The output price and wage are
%   kept explicit throughout; z = p/w is only formed inside the firm-level
%   functions n.m and prof_fn.m.
%
% The input structure contains, among others:
%   p, w    : output price and wage (with p the numeraire, default p = 1)
%   ce      : entry cost, denominated in units of labor (hardcoded)
%   cf_h/l  : operating cost by type, in units of labor
%   svec    : productivity vector
%   F_l,F_h : transition matrices by type
%   G       : startup productivity distribution
%   tol     : convergence criterion
%-----------------------------
% The output is a structure whose fields include (unchanged from the
% original, plus p, w, z and an entry-condition diagnostic):
%   p, w, z            : prices (z = p/w)
%   entry_value_check  : E[V(s0;p,w)] in dollars; equals w*ce at the solved w
%   sstar_ind_l/h      : exit thresholds by type
%   mstar              : mass of (unmeasured) entrants
%   mustar_l/h         : firm-size distribution before exit stage
%   Gnew_l/h, Fnew_l/h : startup dist / transition matrices net of exit
%--------------------------------------------------------------------------
%##########################################################################

function f = stationary(params)

% Get parameters
cf_h    = params.cf_h ;
cf_l    = params.cf_l ;
g       = params.g ;
ns      = params.ns ;
svec    = params.svec ;
F_h     = params.F_h ;
F_l     = params.F_l ;
G       = params.G ;
p       = params.p ;        % output price (numeraire)
w       = params.w ;        % wage (solved from free entry)
ce      = params.ce ;       % entry cost in units of labor (hardcoded input)
mstar   = params.mstar ;
omega_h = params.omega_h;

omega_l = 1 - omega_h;

% Solve for value functions Vh and Vl (in dollars) and thresholds
% sstar_indh and sstar_indl. Values depend on p and w only through z = p/w.
v_h = vfn(F_h, p, w, cf_h, params)  ;
sstar_ind_h = find(v_h) ;
if (isempty(sstar_ind_h))
    results.flag = 0 ;
    f = results ;
    return
else
    sstar_ind_h = sstar_ind_h(1) ;
end

v_l = vfn(F_l, p, w, cf_l, params)  ;
sstar_ind_l = find(v_l) ;
if (isempty(sstar_ind_l))
    results.flag = 0 ;
    f = results ;
    return
else
    sstar_ind_l = sstar_ind_l(1) ;
end

% Diagnostic: expected DOLLAR value of an entrant at the given prices. Free
% entry requires this to equal the dollar entry cost w*ce; solve_wage picks w
% to enforce it (entry_value_check should equal w*ce at the solved wage).
entry_value_check = omega_h*sum(v_h.*G) + omega_l*sum(v_l.*G);

% Degenerate-equilibrium guard: if no firm ever exits (both exit thresholds
% collapse to the first grid point) the prices/ce are outside a sensible
% range. Warn rather than silently return a nonsensical steady state.
if sstar_ind_h == 1 && sstar_ind_l == 1
    warning('stationary:noExit', ...
        ['No firm exits (both exit thresholds at grid index 1) at z = p/w = %.4g. ', ...
         'The equilibrium is likely degenerate; check the hardcoded ce.'], p/w) ;
end

% Find stationary distribution(s) for each mu type

% Start with stationary distribution before exit
Gnew_h = G ;
Gnew_h(1:sstar_ind_h-1) = 0 ;

% New transition matrix to account for exit
Fnew_h = F_h ;
Fnew_h(:,1:sstar_ind_h-1) = zeros(ns, sstar_ind_h-1) ;

% Find m0 such that labor markets clear
I = eye(ns) ;
mustar_h = (I - Fnew_h/(1+g)) \ (omega_h*mstar*Gnew_h) ;
mustar_h(1:sstar_ind_h-1) = 0 ;


% Start with stationary distribution before exit
Gnew_l = G ;
Gnew_l(1:sstar_ind_l-1) = 0 ;

% New transition matrix to account for exit
Fnew_l = F_l ;
Fnew_l(:,1:sstar_ind_l-1) = zeros(ns, sstar_ind_l-1) ;

% Find m0 such that labor markets clear
I = eye(ns) ;
mustar_l = (I - Fnew_l/(1+g)) \ (omega_l*mstar*Gnew_l) ;
mustar_l(1:sstar_ind_l-1) = 0 ;


% calculate employment at each firm

nstar = n(params, svec, p, w);

nstar_h = nstar + cf_h;
nstar_h(1:sstar_ind_h-1) = 0 ;

nstar_l = nstar + cf_l;
nstar_l(1:sstar_ind_l-1) = 0 ;

N = sum(nstar_h.*mustar_h) + sum(nstar_l.*mustar_l) + mstar*ce;

avg_stsize = (omega_h*sum(Gnew_h)*sum(nstar_h.*Gnew_h/sum(Gnew_h)) + omega_l*sum(Gnew_l)*sum(nstar_l.*Gnew_l/sum(Gnew_l)))/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l)) ;
avg_fsize = (sum(nstar_h.*mustar_h)+sum(nstar_l.*mustar_l))/(sum(mustar_h)+sum(mustar_l)) ;
avg_stsizeU = ce + omega_h*sum(nstar_h.*Gnew_h) + omega_l*sum(nstar_l.*Gnew_l); % Denominator of dyn equation: Total employment of 0-year old firms + entry cost.

startup_rate = mstar*(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l))/((sum(mustar_h)+sum(mustar_l))/(1+g)) ;

startup_share_h = omega_h*sum(Gnew_h)/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l));
startup_share_l = omega_l*sum(Gnew_l)/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l));


% -------------------------------------------------------------------------
% Calculate other steady state outcomes
% -------------------------------------------------------------------------

ageub = 250 ;
age=transpose(0:ageub) ;
% Lifecycle Profile for high mu

cohevol_h(:,1) = Gnew_h/sum(Gnew_h) ;
dist_h(:,1) = cohevol_h(:,1)/sum(cohevol_h(:,1)) ;

avg_fsize_age_h = zeros(ageub+1,1); % Preallocate space for speed
avg_fsize_age_h(1) = sum(nstar_h .* dist_h(:,1)) ;
cohevol_emp_h(:,1) = cohevol_h(:,1).*nstar_h ;
for i=2:ageub+1
    cohevol_h(:,i) = Fnew_h * cohevol_h(:,i-1) ;
    cohevol_h(1:sstar_ind_h-1,i) = 0 ;
    dist_h(:,i) = cohevol_h(:,i)/sum(cohevol_h(:,i)) ;
    avg_fsize_age_h(i,1) = sum(nstar_h .* dist_h(:,i)) ;
    cohevol_emp_h(:,i) = cohevol_h(:,i).*nstar_h ;
end
% Unconditional survival probability by age
St_h = sum(cohevol_h(:,1:end)) ;

% Lifecycle Profile for low mu
cohevol_l(:,1) = Gnew_l/sum(Gnew_l) ;
dist_l(:,1) = cohevol_l(:,1)/sum(cohevol_l(:,1)) ;

avg_fsize_age_l = zeros(ageub+1,1); % Preallocate space for speed
avg_fsize_age_l(1) = sum(nstar_l .* dist_l(:,1)) ;
cohevol_emp_l(:,1) = cohevol_l(:,1).*nstar_l ;
for i=2:ageub+1
    cohevol_l(:,i) = Fnew_l * cohevol_l(:,i-1) ;
    cohevol_l(1:sstar_ind_l-1,i) = 0 ;
    dist_l(:,i) = cohevol_l(:,i)/sum(cohevol_l(:,i)) ;
    avg_fsize_age_l(i,1) = sum(nstar_l .* dist_l(:,i)) ;
    cohevol_emp_l(:,i) = cohevol_l(:,i).*nstar_l ;
end

% Unconditional survival probability by age
St_l = sum(cohevol_l(:,1:end)) ;

% Unconditional survival probability by age, for aggregate
St = (startup_share_l * St_l + startup_share_h * St_h);
exit_rate_age =  (St(1:end-1)-St(2:end))./St(1:end-1);
age_weight_l = (startup_share_l*St_l')./(startup_share_l* St_l' + startup_share_h * St_h');
age_weight_h = (startup_share_h*St_h')./(startup_share_l* St_l' + startup_share_h * St_h');
avg_fsize_age = age_weight_l.* avg_fsize_age_l + age_weight_h.* avg_fsize_age_h;


% Conditional growth rate by age
cond_growth_age = avg_fsize_age(2:end)./avg_fsize_age(1:end-1) - 1 ;

% Distribution of firms by age
age_dist = transpose(St)./((1+g).^age) ;

% The last age is interpreted as all firms above that age
age_dist(end) = age_dist(end)/exit_rate_age(end) ;
age_dist = age_dist/sum(age_dist) ;

% Distribution of employment by age
age_empdist = age_dist .* avg_fsize_age ;
age_empdist = age_empdist/sum(age_empdist) ;

% Age-size distribution
szage_mass_h = cohevol_h./repmat(transpose((1+g).^age), [size(cohevol_h,1) 1]) ;
szage_mass_l = cohevol_l./repmat(transpose((1+g).^age), [size(cohevol_l,1) 1]) ;
szageemp_mass_h = szage_mass_h.*repmat(nstar_h, [1 size(szage_mass_h,2)]) ;
szageemp_mass_l = szage_mass_l.*repmat(nstar_l, [1 size(szage_mass_l,2)]) ;

conc_age_h = zeros(ageub+1,1); % Preallocate space for speed
conc_age_l = zeros(ageub+1,1); % Preallocate space for speed
emp_share_age_h = zeros(ageub+1,1); % Preallocate space for speed
emp_share_age_l = zeros(ageub+1,1); % Preallocate space for speed
conc_age = zeros(ageub+1,1); % Preallocate space for speed

% Smallest grid index whose employment reaches 250 (find returns a vector;
% take the first index explicitly to avoid a vector colon-operand warning).
ind_size250_h = find(nstar_h >= 250) ; ind_size250_h = ind_size250_h(1) ;
ind_size250_l = find(nstar_l >= 250) ; ind_size250_l = ind_size250_l(1) ;

for i=1:ageub+1
     conc_age_h(i,1) = sum(szageemp_mass_h(ind_size250_h:end, i)/sum(szageemp_mass_h(:, i))) ;
     conc_age_l(i,1) = sum(szageemp_mass_l(ind_size250_l:end, i)/sum(szageemp_mass_l(:, i))) ;
     emp_share_age_h(i,1) = sum(startup_share_h*szageemp_mass_h(:, i))/(sum(startup_share_h*szageemp_mass_h(:, i))+sum(startup_share_l*szageemp_mass_l(:, i)));
     emp_share_age_l(i,1) = 1 - emp_share_age_h(i,1);
     conc_age(i,1) = conc_age_h(i,1)*emp_share_age_h(i,1) + conc_age_l(i,1)*emp_share_age_l(i,1) ;
end


results.p = p ;
results.w = w ;
results.z = p / w ;                         % relative (real) output price
results.ce = ce ;                           % hardcoded entry cost (input, labor units)
results.entry_value_check = entry_value_check ; % E[V(s0)] in dollars; equals w*ce at solved w
results.N = N ;
results.sstar_ind_h = sstar_ind_h ;
results.sstar_ind_l = sstar_ind_l ;
results.mstar = mstar ;
results.mustar_h = mustar_h ;
results.mustar_l = mustar_l ;
results.nstar = nstar ;
results.omega_h = omega_h;
results.omega_l = omega_l;
results.startup_share_h = startup_share_h;
results.startup_share_l = startup_share_l;
results.nstar_h = nstar_h ;
results.nstar_l = nstar_l ;
results.Gnew_h = Gnew_h ;
results.Gnew_l = Gnew_l ;
results.Fnew_h = Fnew_h ;
results.Fnew_l = Fnew_l ;
results.avg_stsizeU = avg_stsizeU ;
results.avg_stsize = avg_stsize ;
results.avg_fsize = avg_fsize ;
results.startup_rate = startup_rate ;
results.age = age ;
results.age_dist = age_dist ;
results.age_empdist = age_empdist ;
results.exit_rate_age = exit_rate_age' ;
results.avg_fsize_age = avg_fsize_age ;
results.avg_fsize_age_h = avg_fsize_age_h ;
results.avg_fsize_age_l = avg_fsize_age_l ;
results.St = St ;
results.cond_growth_age = cond_growth_age ;
results.conc_age = conc_age ;
results.age_weight_h = age_weight_h;


f = results ;
