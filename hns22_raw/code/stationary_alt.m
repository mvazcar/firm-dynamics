%##########################################################################
%--------------------------------------------------------------------------
% This function solves for the stationary equilibrium for alternative experiments in
% Hopenhayn, Neira and Singhania (2022) and returns the outcome. 
% This code solves for zstar, as opposed to the main code.
% The input structure contains the following variables
% N: number of workers
% ce: entry cost denominated in units of labor
% cf: operating cost denominated in units of labor
% ns: length of productivity vector
% svec: productivity vector
% F_l: transition matrix of low type
% F_h: transition matrix of low type
% G: startup productivity distribution
% tol: convergence criterion
%-----------------------------
% The output is a structure with the following variables
% zstar: price
% sstar_ind_l: exit threshold of low mu types
% sstar_ind_h: exit threshold of high mu types
% mstar: mass of (unmeasured) entrants
% mustar: firm-size distribution before exit stage
% Gnew: startup distribution modified to account for exit
% Fnew_l: transition matrix of low types modified to account for exit
% Fnew_h: transition matrix of high types modified to account for exit
%##########################################################################

function f = stationary_alt(params)

ce = params.ce;
g = params.g ;
mstar = params.mstar ;
ns = params.ns;
svec = params.svec ;
G = params.G ;
F_h = params.F_h ;
F_l = params.F_l ;
omega_h = params.omega_h ;
cf_h = params.cf_h ;
cf_l = params.cf_l ;

% MAIN CHANGE: Take ce as given, solve for zstar that clears entry condition
fun1 = @(zstar)zstar_fun(zstar, params) ;
[zstar, ~] = fzero(fun1, 1.2) ;
params.zstar = zstar ;

% Solve for value functions vh and vl and thresholds sstar_indh and
% sstar_indl
v_h = vfn(F_h, zstar, cf_h, params)  ;
sstar_ind_h = find(v_h) ;
if (isempty(sstar_ind_h))
    stateqbm.flag = 0 ;
    f = results ;
    return
else
    sstar_ind_h = sstar_ind_h(1) ;
end

v_l = vfn(F_l, zstar, cf_l, params)  ;
sstar_ind_l = find(v_l) ;
if (isempty(sstar_ind_l))
    stateqbm.flag = 0 ;
    f = results ;
    return
else
    sstar_ind_l = sstar_ind_l(1) ;
end

omega_l = 1 - omega_h;


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

nstar = n(params, svec, zstar);

nstar_h = nstar + cf_h;
nstar_h(1:sstar_ind_h-1) = 0 ;

nstar_l = nstar + cf_l;
nstar_l(1:sstar_ind_l-1) = 0 ;

N = sum(nstar_h.*mustar_h) + sum(nstar_l.*mustar_l) + mstar*ce;

avg_stsize = (omega_h*sum(Gnew_h)*sum(nstar_h.*Gnew_h/sum(Gnew_h)) + omega_l*sum(Gnew_l)*sum(nstar_l.*Gnew_l/sum(Gnew_l)))/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l)) ;
avg_fsize = (sum(nstar_h.*mustar_h)+sum(nstar_l.*mustar_l))/(sum(mustar_h)+sum(mustar_l)) ;

startup_rate = mstar*(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l))/((sum(mustar_h)+sum(mustar_l))/(1+g)) ;

% 250+ cutoff
ind_size250_h = find(nstar_h >= 250) ;
ind_size250_l = find(nstar_l >= 250) ;

% Concentration by type
conc_h = sum(mustar_h(ind_size250_h:end).*nstar_h(ind_size250_h:end))/sum(mustar_h.*nstar_h);
conc_l = sum(mustar_l(ind_size250_l:end).*nstar_l(ind_size250_l:end))/sum(mustar_l.*nstar_l);

emp_share_h = sum(mustar_h.*nstar_h)/(sum(mustar_h.*nstar_h)+sum(mustar_l.*nstar_l));
emp_share_l = sum(mustar_l.*nstar_l)/(sum(mustar_h.*nstar_h)+sum(mustar_l.*nstar_l));

% Aggregate Concentration 
conc = conc_h*emp_share_h + conc_l*emp_share_l;

startup_rate_h = mstar*omega_h*(sum(Gnew_h))/(sum(mustar_h/(1+g))) ; % Number of new firms is mstar*omegah*(sum(Gnewh))
startup_rate_l = mstar*omega_l*(sum(Gnew_l))/(sum(mustar_l/(1+g)));

startup_share_h = omega_h*sum(Gnew_h)/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l));
startup_share_l = omega_l*sum(Gnew_l)/(omega_h*sum(Gnew_h)+omega_l*sum(Gnew_l));

% -------------------------------------------------------------------------
% Calculate other steady state outcomes: 
% 1. firm size by age 2. exit rate by age 3. aging 4. concentration and 
% 5. distributions
% -------------------------------------------------------------------------

ageub = 250 ;
age=transpose(0:ageub) ;


% Lifecycle Profile for high mu

cohevol_h(:,1) = Gnew_h/sum(Gnew_h) ;
dist_h(:,1) = cohevol_h(:,1)/sum(cohevol_h(:,1)) ;

avg_fsize_age_h(1) = sum(nstar_h .* dist_h(:,1)) ;
cohevol_emp_h(:,1) = cohevol_h(:,1).*nstar_h ;
for i=2:ageub+1
        
    cohevol_h(:,i) = Fnew_h * cohevol_h(:,i-1) ;
    cohevol_h(1:sstar_ind_h-1,i) = 0 ;
    dist_h(:,i) = cohevol_h(:,i)/sum(cohevol_h(:,i)) ;

    exit_rate_age_h(i-1,1) = (sum(cohevol_h(:,i-1)) - sum(cohevol_h(:,i)))/sum(cohevol_h(:,i-1)) ;
    exit_rate_age_contemp_h(i-1,1) = (sum(cohevol_h(:,i-1)) - sum(cohevol_h(:,i)))/sum(cohevol_h(:,i)) ;
    avg_fsize_age_h(i,1) = sum(nstar_h .* dist_h(:,i)) ;

    cohevol_emp_h(:,i) = cohevol_h(:,i).*nstar_h ;
end
% Unconditional survival probability by age
St_h = sum(cohevol_h(:,1:end)) ;

% Lifecycle Profile for low mu

cohevol_l(:,1) = Gnew_l/sum(Gnew_l) ;
dist_l(:,1) = cohevol_l(:,1)/sum(cohevol_l(:,1)) ;

avg_fsize_age_l(1) = sum(nstar_l .* dist_l(:,1)) ;
cohevol_emp_l(:,1) = cohevol_l(:,1).*nstar_l ;
for i=2:ageub+1
       
    cohevol_l(:,i) = Fnew_l * cohevol_l(:,i-1) ;
    cohevol_l(1:sstar_ind_l-1,i) = 0 ;
    dist_l(:,i) = cohevol_l(:,i)/sum(cohevol_l(:,i)) ;

    exit_rate_age_l(i-1,1) = (sum(cohevol_l(:,i-1)) - sum(cohevol_l(:,i)))/sum(cohevol_l(:,i-1)) ;
    exit_rate_age_contemp_l(i-1,1) = (sum(cohevol_l(:,i-1)) - sum(cohevol_l(:,i)))/sum(cohevol_l(:,i)) ;
    avg_fsize_age_l(i,1) = sum(nstar_l .* dist_l(:,i)) ;
    cohevol_emp_l(:,i) = cohevol_l(:,i).*nstar_l ;
end

% Unconditional survival probability by age
St_l = sum(cohevol_l(:,1:end)) ;

% Unconditional survival probability by age, for aggregate
St = (startup_share_l * St_l + startup_share_h * St_h);
exit_rate_age =  (St(1:end-1)-St(2:end))./St(1:end-1);
exit_rate_age_contemp =  (St(1:end-1)-St(2:end))./St(2:end);
age_weight_l = (startup_share_l*St_l')./(startup_share_l* St_l' + startup_share_h * St_h');
age_weight_h = (startup_share_h*St_h')./(startup_share_l* St_l' + startup_share_h * St_h');
avg_fsize_age = age_weight_l.* avg_fsize_age_l + age_weight_h.* avg_fsize_age_h;

% Conditional growth rate by age
cond_growth_age = avg_fsize_age(2:end)./avg_fsize_age(1:end-1) - 1 ;

% Firm size distribution in stat eqbm
fsize_dist_h = mustar_h/sum(mustar_h) ;
fsize_dist_l = mustar_l/sum(mustar_l) ;

% Firm employment distribution in stat eqbm

femp_dist_h = mustar_h.*nstar_h;
femp_dist_h = femp_dist_h/sum(femp_dist_h) ;

femp_dist_l = mustar_l.*nstar_l;
femp_dist_l = femp_dist_l/sum(femp_dist_l) ;


% Distribution of firms by age 
age_dist = transpose(St)./((1+g).^age) ;
% The last age is interpreted as all firms above that age
age_dist(end) = age_dist(end)/exit_rate_age(end) ;
age_dist = age_dist/sum(age_dist) ;
% Calculate fraction of 11+ firms in stationary eqbm
ind_age11 = find(age == 11) ;
frac11plus = sum(age_dist(ind_age11:end)) ;

% Calculate average age in the economy, to show aging
avg_age = sum(age .* age_dist) ;

% Calculate average age of Above 25 group, to see how group is aging
avg_age_above25 = sum(age(27:end) .* age_dist(27:end))/sum(age_dist(27:end)) ;


% Distribution of employment by age
age_empdist = age_dist .* avg_fsize_age ;
age_empdist = age_empdist/sum(age_empdist) ;
frac11plus_emp = sum(age_empdist(ind_age11:end)) ;

% Distribution of employment by age by type
age_dist_h = transpose(St_h)./((1+g).^age) ;
age_dist_h(end) = age_dist_h(end)/exit_rate_age_h(end) ;
age_dist_h = age_dist_h/sum(age_dist_h) ;
age_empdist_h = age_dist_h .* avg_fsize_age_h ;
age_empdist_h = age_empdist_h/sum(age_empdist_h) ;
frac11plus_emp_h = sum(age_empdist_h(ind_age11:end)) ;

age_dist_l = transpose(St_l)./((1+g).^age) ;
age_dist_l(end) = age_dist_l(end)/exit_rate_age_h(end) ;
age_dist_l = age_dist_l/sum(age_dist_l) ;
age_empdist_l = age_dist_l .* avg_fsize_age_h ;
age_empdist_l = age_empdist_l/sum(age_empdist_l) ;
frac11plus_emp_l = sum(age_empdist_l(ind_age11:end)); 

% Age-size distribution
szage_mass_h = cohevol_h./repmat(transpose((1+g).^age), [size(cohevol_h,1) 1]) ;
szageemp_mass_h = szage_mass_h.*repmat(nstar_h, [1 size(szage_mass_h,2)]) ;
szage_dist_h = szage_mass_h/sum(sum(szage_mass_h)) ;
szageemp_dist_h = szageemp_mass_h/sum(sum(szageemp_mass_h)) ;

szage_mass_l = cohevol_l./repmat(transpose((1+g).^age), [size(cohevol_l,1) 1]) ;
szageemp_mass_l = szage_mass_l.*repmat(nstar_l, [1 size(szage_mass_l,2)]) ;
szage_dist_l = szage_mass_l/sum(sum(szage_mass_l)) ;
szageemp_dist_l = szageemp_mass_l/sum(sum(szageemp_mass_l)) ;

szage_mass = startup_share_h*szage_mass_h + startup_share_l*szage_mass_l;
szageemp_mass = szage_mass.*repmat(nstar, [1 size(szage_mass,2)]) ;
szage_dist = szage_mass/sum(sum(szage_mass)) ;
szageemp_dist = szageemp_mass/sum(sum(szageemp_mass)) ;


for i=1:ageub+1
     conc_age_h(i,1) = sum(szageemp_mass_h(ind_size250_h:end, i)/sum(szageemp_mass_h(:, i))) ;
     conc_age_l(i,1) = sum(szageemp_mass_l(ind_size250_l:end, i)/sum(szageemp_mass_l(:, i))) ;
     emp_share_age_h(i,1) = sum(startup_share_h*szageemp_mass_h(:, i))/(sum(startup_share_h*szageemp_mass_h(:, i))+sum(startup_share_l*szageemp_mass_l(:, i)));
     emp_share_age_l(i,1) = 1 - emp_share_age_h(i,1);
     conc_age(i,1) = conc_age_h(i,1)*emp_share_age_h(i,1) + conc_age_l(i,1)*emp_share_age_l(i,1) ;
end


% Verify agg. concentration adds up
conc_fromAge = sum(conc_age.*age_empdist);

% Verify agg. firm size adds up
afs_fromAge = sum(avg_fsize_age.*age_dist); 

% Verify agg. exit adds up
exit_fromAge = sum(exit_rate_age'.*age_dist(1:end-1)/sum(age_dist(1:end-1)));



%% Verify that startup rate equals labor force growth rate plus exit rate
% We will do so by feeding a time series for labor that grows at a
% constant rate and solving for a fill in the blanks equilibrium

tau = 100 ;
t = transpose(0:tau) ;
Ns(1) = N ;
for i=2:length(t)
    Ns(i,1) = Ns(i-1)*(1+g) ;
end


% Specify the matrix to store firm-size distribution for each t, high type
MU_h = zeros(length(mustar_h), length(t)) ;
MU_h(:,1) = mustar_h ;

MUemp_h = zeros(length(mustar_h), length(t)) ;
MUemp_h(:,1) = MU_h(:,1).* nstar_h ;

% Ith = zeros(length(t), length(mustar)) ; % incumbent masses
It_h(:,1) = Fnew_h*(mustar_h/(1+g)) ;
It_h(1:sstar_ind_h-1,1) = 0 ;
mtUn_h(1) = omega_h*mstar ;
MtUn_h(1) = sum(It_h(:,1)) + mtUn_h(1) ;
Mt_h(1) = sum(MU_h(:,1)) ;
Xt_h(1) = Mt_h(1)/(1+g) - sum(It_h(:,1)) ;
% Compute the evolution of the age distribution
MU_age_h = zeros(length(age), length(t)) ;
MU_age_h(:,1) = age_dist_h*Mt_h(1) ;
frac11plus_h(1) = sum(MU_age_h(ind_age11:end,1))/sum(MU_age_h(:,1)) ; 
% Employment distribution by age
MUemp_age_h = zeros(length(age), length(t)) ;
MUemp_age_h(:,1) = age_empdist_h*Ns(1) ;
frac11plus_emp_h(1) = sum(MUemp_age_h(ind_age11:end,1))/sum(MUemp_age_h(:,1)) ; 


% Specify the matrix to store firm-size distribution for each t, low type
MU_l = zeros(length(mustar_l), length(t)) ;
MU_l(:,1) = mustar_l ;

MUemp_l = zeros(length(mustar_l), length(t)) ;
MUemp_l(:,1) = MU_l(:,1).* nstar_l;

Mt_l(1) = sum(MU_l(:,1)) ;

It_l(:,1) = Fnew_l*(mustar_l/(1+g)) ;
It_l(1:sstar_ind_l-1,1) = 0 ;

mtUn_l(1) = omega_l*mstar ;
MtUn_l(1) = sum(It_l(:,1)) + mtUn_l(1) ;

Xt_l(1) = Mt_l(1)/(1+g) - sum(It_l(:,1)) ;
% Compute the evolution of the age distribution
MU_age_l = zeros(length(age), length(t)) ;
MU_age_l(:,1) = age_dist_h*Mt_l(1) ;
frac11plus_l(1) = sum(MU_age_l(ind_age11:end,1))/sum(MU_age_l(:,1)) ; 
% Employment distribution by age
MUemp_age_l = zeros(length(age), length(t)) ;
MUemp_age_l(:,1) = age_empdist_l*Ns(1) ;
frac11plus_emp_l(1) = sum(MUemp_age_l(ind_age11:end,1))/sum(MUemp_age_l(:,1)) ; 

% Some aggregates
Mt(1) = Mt_l(1)+Mt_h(1);
AFSt(1) = (sum(nstar_l .* MU_l(:,1))+(sum(nstar_h .* MU_h(:,1))))/(Mt_l(1)+Mt_h(1)) ;
startup_rate_t(1) = startup_rate ;
startup_rate_contemp_t(1) = startup_rate ;
exit_rate_t(1) = (Xt_l(1)+Xt_h(1))/(Mt_l(1)+Mt_h(1))*(1+g) ;
exit_rate_contemp_t(1) = (Xt_l(1)+Xt_h(1))/(Mt_l(1)+Mt_h(1)) ;

avg_stsizeU = ce + omega_h*sum(nstar_h.*Gnew_h) + omega_l*sum(nstar_l.*Gnew_l); % Denominator of dyn equation: Total employment of 0-year old firms + entry cost.

% % Now let's iterate forward by time-periods and calculate evolution of M
for i=2:length(t)
   It_h(:,i) = Fnew_h*MU_h(:,i-1) ;
   It_h(1:sstar_ind_h-1,i) = 0;

   It_l(:,i) = Fnew_l*MU_l(:,i-1) ;
   It_l(1:sstar_ind_l-1,i) = 0;

   It(:,i) = It_l(:,i)+ It_h(:,i); 
   
   Nd_inc = sum(nstar_h .* It_h(:,i)) + sum(nstar_l .* It_l(:,i));
   Nd_res = Ns(i) - Nd_inc ;
   
   m = Nd_res/avg_stsizeU;
   
   % Calculate evolution for high types
   mtUn_h(i,1) = omega_h*m ;  
   mth(i) = mtUn_h(i)*sum(Gnew_h) ;
   MU_h(:,i) = mtUn_h(i,1)*Gnew_h +  It_h(:,i) ;
   MUemp_h(:,i) = MU_h(:,i).* nstar_h;
   Mt_h(i,1) = sum(MU_h(:,i)) ;
   Xt_h(i) = Mt_h(i-1) - sum(It_h(:,i)) ;
 
   % Calculate evolution for low types
   mtUn_l(i,1) = omega_l*m ;
   mtl(i) = mtUn_l(i)*sum(Gnew_l) ;
   MU_l(:,i) = mtUn_l(i,1)*Gnew_l +  It_l(:,i) ;
   MUemp_l(:,i) = MU_l(:,i).* nstar_l;
   Mt_l(i,1) = sum(MU_l(:,i)) ;
   Xt_l(i) = Mt_l(i-1) - sum(It_l(:,i)) ;
   
   Mt(i) = Mt_l(i)+Mt_h(i);
   mt(i) = mtl(i)+mth(i);
   Xt(i) = Xt_l(i)+Xt_h(i);
   
 
   % measured aggregates
    AFSt(i,1) = (sum(nstar_l.* MU_l(:,i))+(sum(nstar_h.* MU_h(:,i))))/Mt(i) ;
    startup_rate_t(i) = mt(i)/Mt(i-1);
    exit_rate_t(i) = Xt(i)/Mt(i-1) ;
   
end


results.startup_rate = startup_rate;
results.exit_rate = exit_rate_t(1);
results.avg_fsize = avg_fsize;
results.conc = conc;

results.exit_rate_age = exit_rate_age;
results.avg_fsize_age = avg_fsize_age;
results.conc_age = conc_age;

results.age_dist = age_dist;
results.age_empdist = age_empdist;

f = results ;
