%--------------------------------------------------------------------------
% This function takes in the estimated model and feeds labor force growth
% rates observed in the data through the model. It returns a structure
% containing various statistics.
% Inputs:
%   filename: contains the full path to the csv file containing the data.
%   The data are stored in a uniform format: year, civilian labor force
%--------------------------------------------------------------------------
function f = feed_data(params, eqbm, filename)

% Read data on labor force growth
lf_data = csvread(filename, 1, 0) ;
year = lf_data(1:end,1) ;  % Read year
year = [year(1) - 1; year] ; % Add an extra year in the beginning
lf_growth = lf_data(1:end,2) ; % Read lf growth data
lf_growth = [1; lf_growth] ; % Normalize lf growth in initial year to 1+g
Ns = eqbm.N*cumprod(lf_growth) ; % Calculate time-series of labor force

% Number of periods to simulate forward
t = 1:length(year) ;

%------------------------------------------------------------------
% Firm size dist: low type
%------------------------------------------------------------------
% Size distribution
MU_l = zeros(length(t), length(eqbm.mustar_l)) ; % firm masses
MU_l(1,:) = eqbm.mustar_l ; % initialize

% Employment distribution
MUemp_l = zeros(length(t), length(eqbm.mustar_l)) ; % Total workers by size
MUemp_l(1,:) = MU_l(1,:).*transpose(eqbm.nstar_l) ; % initialize

It_l = zeros(length(t), length(eqbm.mustar_l)) ; % incumbent masses
It_l(1,:) = eqbm.Fnew_l*(eqbm.mustar_l/(1+params.g)) ; % initialize
It_l(1,1:eqbm.sstar_ind_l-1) = 0 ; % get rid of firms below exit threshold

%------------------------------------------------------------------
% Firm size dist: high type
%------------------------------------------------------------------
% Size distribution
MU_h = zeros(length(t), length(eqbm.mustar_h)) ; % firm masses
MU_h(1,:) = eqbm.mustar_h ; % initialize

% Employment distribution
MUemp_h = zeros(length(t), length(eqbm.mustar_h)) ; % Total workers by size
MUemp_h(1,:) = MU_h(1,:).*transpose(eqbm.nstar_h) ; % initialize

It_h = zeros(length(t), length(eqbm.mustar_h)) ; % incumbent masses
It_h(1,:) = eqbm.Fnew_h*(eqbm.mustar_h/(1+params.g)) ; % initialize
It_h(1,1:eqbm.sstar_ind_h-1) = 0 ; % get rid of firms below exit threshold

%------------------------------------------------------------------
% Firm size dist: Aggregates
%------------------------------------------------------------------

% Fraction of firms with 250+ employees in masses and by employment
empthres = 250 ; % threshold that defines large firms
frac250plus = NaN(length(t),1) ; 
frac250plus_emp = NaN(length(t),1) ;
frac1kplus_emp = NaN(length(t),1) ;
frac10kplus_emp = NaN(length(t),1) ;

frac250plus(1) = (sum(MU_h(1,eqbm.nstar_h >= empthres))+ sum(MU_l(1,eqbm.nstar_l >= empthres)))/(sum(MU_l(1,:))+sum(MU_h(1,:)));

frac250plus_emp(1) = (sum(MUemp_h(1,eqbm.nstar_h >= empthres))+ sum(MUemp_l(1,eqbm.nstar_l >= empthres)))/(sum(MUemp_l(1,:))+sum(MUemp_h(1,:)));

frac1kplus_emp(1) = (sum(MUemp_h(1,eqbm.nstar_h >= 999))+ sum(MUemp_l(1,eqbm.nstar_l >= 999)))/(sum(MUemp_l(1,:))+sum(MUemp_h(1,:)));

frac10kplus_emp(1) = (sum(MUemp_h(1,eqbm.nstar_h >= 9999))+ sum(MUemp_l(1,eqbm.nstar_l >= 9999)))/(sum(MUemp_l(1,:))+sum(MUemp_h(1,:)));

% Fraction 1 to 4
frac1to4 = NaN(length(t),1) ; 
frac1to4(1) = (sum(MU_h(1,eqbm.nstar_h <= 4)) + sum(MU_l(1,eqbm.nstar_l <= 4)))/(sum(MU_l(1,:))+sum(MU_h(1,:)));    

%------------------------------------------------------------------
% Unmeasured entry-exit aggregates
%------------------------------------------------------------------
mtUn_l = NaN(length(t),1) ; % Unmeasured mass of entrants
mtUn_l(1) = eqbm.omega_l*params.mstar ; % initialize


mtUn_h = NaN(length(t),1) ; % Unmeasured mass of entrants
mtUn_h(1) = eqbm.omega_h*params.mstar ; % initialize

mtUn = NaN(length(t),1) ; % Unmeasured mass of entrants
mtUn(1) = params.mstar ; % initialize
MtUn = NaN(length(t),1) ; % Mass of firms
MtUn(1) = sum(It_l(1,:) + It_h(1,:)) + mtUn(1) ; % initialize

%------------------------------------------------------------------
% Measured entry-exit aggregates
%------------------------------------------------------------------
mt_l = NaN(length(t),1) ; % Mass of entrants
mt_l(1) = eqbm.omega_l*params.mstar*sum(eqbm.Gnew_l) ; % initialize
Mt_l = NaN(length(t),1) ; % Mass of firms
Mt_l(1) = sum(MU_l(1,:)) ; % initialize
Xt_l = NaN(length(t),1) ; % Exit mass
Xt_l(1) = Mt_l(1)/(1+params.g) - sum(It_l(1,:)) ; % initialize 


mt_h = NaN(length(t),1) ; % Mass of entrants
mt_h(1) = eqbm.omega_h*params.mstar*sum(eqbm.Gnew_h) ; % initialize
Mt_h = NaN(length(t),1) ; % Mass of firms
Mt_h(1) = sum(MU_h(1,:)) ; % initialize
Xt_h = NaN(length(t),1) ; % Exit mass
Xt_h(1) = Mt_h(1)/(1+params.g) - sum(It_h(1,:)) ; % initialize 


% Now calculate some aggregates from low and high types
Mt = NaN(length(t),1);
Mt(1) = Mt_l(1)+Mt_h(1);

E_l = NaN(length(t),1);
E_l(1) = sum(transpose(eqbm.nstar_l) .* MU_l(1,:)); % Total Employment by small mu firms

E_h = NaN(length(t),1);
E_h(1) = sum(transpose(eqbm.nstar_h) .* MU_h(1,:)); % Total Employment by large mu firms

E = NaN(length(t),1);
E(1) = E_l(1) + E_h(1); % Total Employment

mt = NaN(length(t),1);
mt(1) = mt_l(1)+mt_h(1); % Total number of measured startups

Xt = NaN(length(t),1);
Xt(1) = Xt_l(1)+Xt_h(1); % Total Exits

mtUn = NaN(length(t),1);
mtUn(1) = eqbm.mstar; % Total number of measured startups

AFSt = NaN(length(t),1) ; % Average firm size
AFSt(1) = E(1)/Mt(1) ; % initialize
startup_rate = NaN(length(t),1) ; % Startup rate
startup_rate(1) = mt(1)/(Mt(1)/(1 + params.g)) ; % initialize
exit_rate = NaN(length(t),1) ; % Exit rate
exit_rate(1) = Xt(1)/(Mt(1)/(1 + params.g)) ; % initialize


%------------------------------------------------------------------
% Firm age
%------------------------------------------------------------------
MU_age = zeros(length(t), length(eqbm.age)) ; % firm masses by age
MU_age(1,:) = transpose(eqbm.age_dist)*Mt(1) ;

MUemp_age = zeros(length(t), length(eqbm.age)) ; % number of workers by firm age
MUemp_age(1,:) = transpose(eqbm.age_empdist)*(sum(MUemp_h(1,:))+ sum(MUemp_l(1,:))) ;

ind_age6 = find(eqbm.age == 6) ;
ind_age10 = find(eqbm.age == 10) ;
ind_age11 = find(eqbm.age == 11) ;
ind_age15 = find(eqbm.age == 15) ;
ind_age16 = find(eqbm.age == 16) ;
ind_age20 = find(eqbm.age == 20) ;
ind_age21 = find(eqbm.age == 21) ;
ind_age25 = find(eqbm.age == 25) ;
ind_age26 = find(eqbm.age == 26) ;

avg_fsize_age_6_10 = NaN(length(t),1) ; 
avg_fsize_age_11_15 = NaN(length(t),1) ; 
avg_fsize_age_16_20 = NaN(length(t),1) ; 
avg_fsize_age_21_25 = NaN(length(t),1) ; 
avg_fsize_age_above25 = NaN(length(t),1) ; 

avg_fsize_age_6_10(1) = mean(sum(eqbm.avg_fsize_age(ind_age6:ind_age10).*MU_age(1,ind_age6:ind_age10)/sum(MU_age(1,ind_age6:ind_age10)))); 
avg_fsize_age_11_15(1) = mean(sum(eqbm.avg_fsize_age(ind_age11:ind_age15).*MU_age(1,ind_age11:ind_age15)/sum(MU_age(1,ind_age11:ind_age15)))); 
avg_fsize_age_16_20(1) = mean(sum(eqbm.avg_fsize_age(ind_age16:ind_age20).*MU_age(1,ind_age16:ind_age20)/sum(MU_age(1,ind_age16:ind_age20)))); 
avg_fsize_age_21_25(1) = mean(sum(eqbm.avg_fsize_age(ind_age21:ind_age25).*MU_age(1,ind_age21:ind_age25)/sum(MU_age(1,ind_age21:ind_age25)))); 
avg_fsize_age_above25(1) = mean(sum(eqbm.avg_fsize_age(ind_age26:end)'.*MU_age(1,ind_age26:end)/sum(MU_age(1,ind_age26:end)))); 

exit_rate_age_6_10 = NaN(length(t),1) ; 
exit_rate_age_11_15 = NaN(length(t),1) ; 
exit_rate_age_16_20 = NaN(length(t),1) ; 
exit_rate_age_21_25 = NaN(length(t),1) ; 
exit_rate_age_above25 = NaN(length(t),1) ; 

exit_rate_age_6_10(1) = mean(sum(eqbm.exit_rate_age(ind_age6-1:ind_age10-1)'.*MU_age(1,ind_age6-1:ind_age10-1)/sum(MU_age(1,ind_age6-1:ind_age10-1))));
exit_rate_age_11_15(1) = mean(sum(eqbm.exit_rate_age(ind_age11-1:ind_age15-1)'.*MU_age(1,ind_age11-1:ind_age15-1)/sum(MU_age(1,ind_age11-1:ind_age15-1))));
exit_rate_age_16_20(1) = mean(sum(eqbm.exit_rate_age(ind_age16-1:ind_age20-1)'.*MU_age(1,ind_age16-1:ind_age20-1)/sum(MU_age(1,ind_age16-1:ind_age20-1))));
exit_rate_age_21_25(1) = mean(sum(eqbm.exit_rate_age(ind_age21-1:ind_age25-1)'.*MU_age(1,ind_age21-1:ind_age25-1)/sum(MU_age(1,ind_age21-1:ind_age25-1))));
exit_rate_age_above25(1) = mean(sum(eqbm.exit_rate_age(ind_age25:end)'.*MU_age(1,ind_age25:end-1)/sum(MU_age(1,ind_age25:end-1))));

conc_age_6_10 = NaN(length(t),1) ; 
conc_age_11_15 = NaN(length(t),1) ; 
conc_age_16_20 = NaN(length(t),1) ; 
conc_age_21_25 = NaN(length(t),1) ; 
conc_age_above25 = NaN(length(t),1) ; 

conc_age_6_10(1) = mean(sum(eqbm.conc_age(ind_age6:ind_age10)'.*MUemp_age(1,ind_age6:ind_age10)/sum(MUemp_age(1,ind_age6:ind_age10))));
conc_age_11_15(1) = mean(sum(eqbm.conc_age(ind_age11:ind_age15)'.*MUemp_age(1,ind_age11:ind_age15)/sum(MUemp_age(1,ind_age11:ind_age15))));
conc_age_16_20(1) = mean(sum(eqbm.conc_age(ind_age16:ind_age20)'.*MUemp_age(1,ind_age16:ind_age20)/sum(MUemp_age(1,ind_age16:ind_age20))));
conc_age_21_25(1) = mean(sum(eqbm.conc_age(ind_age21:ind_age25)'.*MUemp_age(1,ind_age21:ind_age25)/sum(MUemp_age(1,ind_age21:ind_age25))));
conc_age_above25(1) = mean(sum(eqbm.conc_age(ind_age26:end)'.*MUemp_age(1,ind_age26:end)/sum(MUemp_age(1,ind_age26:end))));


%------------------------------------------------------------------
% Labor share
%------------------------------------------------------------------
% revenue for each firm size
rev_firms = (params.zstar*exp(transpose(params.svec))*params.alpha^params.alpha).^(1/(1-params.alpha)) ;
aggrev = NaN(length(t),1) ; % aggregate revenue
aggrev(1) = sum(MU_h(1,:).*rev_firms + MU_l(1,:).*rev_firms) ;

vw_firms_l = zeros(length(t), length(params.svec)) ; % value added weights
vw_firms_l(1,:) = MU_l(1,:).*rev_firms/aggrev(1) ;

vw_firms_h = zeros(length(t), length(params.svec)) ; % value added weights
vw_firms_h(1,:) = MU_h(1,:).*rev_firms/aggrev(1) ;

lshare_firms_l = transpose(eqbm.nstar_l)./rev_firms ; % labor share by firm size
lshare_firms_h = transpose(eqbm.nstar_h)./rev_firms ; % labor share by firm size

agglshare = NaN(length(t), 1) ; % aggregate labor share
agglshare(1) = sum(lshare_firms_h.*vw_firms_h(1,:))+ sum(lshare_firms_l.*vw_firms_l(1,:)) ;


%------------------------------------------------------------------
% Simulate economy forward
%------------------------------------------------------------------
for i=2:length(t)
    % Fill-in-the-blanks equilibrium
    It_l(i,:) = eqbm.Fnew_l*transpose(MU_l(i-1,:)) ; % measure of incumbents
    It_l(i,1:eqbm.sstar_ind_l-1) = 0 ; % account for incumbent exit
    
    It_h(i,:) = eqbm.Fnew_h*transpose(MU_h(i-1,:)) ; % measure of incumbents
    It_h(i,1:eqbm.sstar_ind_h-1) = 0 ; % account for incumbent exit
    
    Nd_inc = sum(transpose(eqbm.nstar_h) .* It_h(i,:)) + sum(transpose(eqbm.nstar_l) .* It_l(i,:));
    Nd_res = Ns(i) - Nd_inc ; % residual labor supply

   
    % unmeasured variables
    mtUn(i) = Nd_res/eqbm.avg_stsizeU ; % startups
    MtUn(i) = sum(It_l(i,:) + It_h(i,:)) + mtUn(i) ; % firms
   
    % firm size
    mtUn_l(i) = eqbm.omega_l*mtUn(i) ; % startups
    mtUn_h(i) = eqbm.omega_h*mtUn(i) ; % startups
    
    MU_l(i,:) = mtUn_l(i,1)*transpose(eqbm.Gnew_l) +  It_l(i,:) ;    
    MU_h(i,:) = mtUn_h(i,1)*transpose(eqbm.Gnew_h) +  It_h(i,:) ;
    
    MUemp_l(i,:) = MU_l(i,:).* transpose(eqbm.nstar_l) ;
    MUemp_h(i,:) = MU_h(i,:).* transpose(eqbm.nstar_h) ;
    
    frac250plus(i,1) = (sum(MU_h(i,eqbm.nstar_h >= empthres))+ sum(MU_l(i,eqbm.nstar_l >= empthres)))/(sum(MU_l(i,:))+sum(MU_h(i,:)));
    frac250plus_emp(i,1) = (sum(MUemp_h(i,eqbm.nstar_h >= empthres))+ sum(MUemp_l(i,eqbm.nstar_l >= empthres)))/(sum(MUemp_l(i,:))+sum(MUemp_h(i,:)));
    frac1kplus_emp(i,1) = (sum(MUemp_h(i,eqbm.nstar_h >= 999))+ sum(MUemp_l(i,eqbm.nstar_l >= 999)))/(sum(MUemp_l(i,:))+sum(MUemp_h(i,:)));
    frac10kplus_emp(i,1) = (sum(MUemp_h(i,eqbm.nstar_h >= 9999))+ sum(MUemp_l(i,eqbm.nstar_l >= 9999)))/(sum(MUemp_l(i,:))+sum(MUemp_h(i,:)));

    frac1to4(i,1) = (sum(MU_h(i,eqbm.nstar_h <= 4)) + sum(MU_l(i,eqbm.nstar_l <= 4)))/(sum(MU_l(i,:))+sum(MU_h(i,:)));    

    % measured aggregates 
    mt_l(i) = eqbm.omega_l*mtUn(i)*sum(eqbm.Gnew_l) ; % 
    mt_h(i) = eqbm.omega_h*mtUn(i)*sum(eqbm.Gnew_h) ; % 
    mt(i) = mt_l(i) + mt_h(i);
    
    Mt_l(i) = sum(MU_l(i,:)) ; % Total number of low firms
    Mt_h(i) = sum(MU_h(i,:)) ; % Total number of high firms
    Mt(i) = Mt_l(i)+Mt_h(i);

    E_l(i) = sum(transpose(eqbm.nstar_l) .* MU_l(i,:)); % Total Employment by small mu firms
    E_h(i) = sum(transpose(eqbm.nstar_h) .* MU_h(i,:)); % Total Employment by large mu firms
    E(i) = E_l(i) + E_h(i); % Total Employment
    
    Xt(i) = Mt(i-1) - sum(It_l(i,:) + It_h(i,:)) ;
    AFSt(i) =(sum(transpose(eqbm.nstar_l).*MU_l(i,:))+sum(transpose(eqbm.nstar_h).*MU_h(i,:)))/Mt(i) ;
    startup_rate(i) = mt(i)/Mt(i-1) ;
    exit_rate(i) = Xt(i)/Mt(i-1) ;
    
    % Firm age
    % Assuming that exit rate stays constant after oldest age
    Inc_age = MU_age(i-1, 1:end-1).*transpose(1 - eqbm.exit_rate_age) ;
    Inc_age(end) = Inc_age(end) + MU_age(i-1,end)*(1 - eqbm.exit_rate_age(end)) ;
    MU_age(i,:) = [mt(i), Inc_age] ;
    MUemp_age(i,:) = MU_age(i,:).*transpose(eqbm.avg_fsize_age) ;
    
    avg_fsize_age_6_10(i) = mean(sum(eqbm.avg_fsize_age(ind_age6:ind_age10).*MU_age(i,ind_age6:ind_age10)/sum(MU_age(i,ind_age6:ind_age10)))); 
    avg_fsize_age_11_15(i) = mean(sum(eqbm.avg_fsize_age(ind_age11:ind_age15).*MU_age(i,ind_age11:ind_age15)/sum(MU_age(i,ind_age11:ind_age15)))); 
    avg_fsize_age_16_20(i) = mean(sum(eqbm.avg_fsize_age(ind_age16:ind_age20).*MU_age(i,ind_age16:ind_age20)/sum(MU_age(i,ind_age16:ind_age20)))); 
    avg_fsize_age_21_25(i) = mean(sum(eqbm.avg_fsize_age(ind_age21:ind_age25).*MU_age(i,ind_age21:ind_age25)/sum(MU_age(i,ind_age21:ind_age25)))); 
    avg_fsize_age_above25(i) = mean(sum(eqbm.avg_fsize_age(ind_age26:end)'.*MU_age(i,ind_age26:end)/sum(MU_age(i,ind_age26:end)))); 
    
    exit_rate_age_6_10(i) = mean(sum(eqbm.exit_rate_age(ind_age6-1:ind_age10-1)'.*MU_age(i,ind_age6-1:ind_age10-1)/sum(MU_age(i,ind_age6-1:ind_age10-1))));
    exit_rate_age_11_15(i) = mean(sum(eqbm.exit_rate_age(ind_age11-1:ind_age15-1)'.*MU_age(i,ind_age11-1:ind_age15-1)/sum(MU_age(i,ind_age11-1:ind_age15-1))));
    exit_rate_age_16_20(i) = mean(sum(eqbm.exit_rate_age(ind_age16-1:ind_age20-1)'.*MU_age(i,ind_age16-1:ind_age20-1)/sum(MU_age(i,ind_age16-1:ind_age20-1))));
    exit_rate_age_21_25(i) = mean(sum(eqbm.exit_rate_age(ind_age21-1:ind_age25-1)'.*MU_age(i,ind_age21-1:ind_age25-1)/sum(MU_age(i,ind_age21-1:ind_age25-1))));
    exit_rate_age_above25(i) = mean(sum(eqbm.exit_rate_age(ind_age25:end)'.*MU_age(i,ind_age25:end-1)/sum(MU_age(i,ind_age25:end-1))));
    
    conc_age_6_10(i) = mean(sum(eqbm.conc_age(ind_age6:ind_age10)'.*MUemp_age(i,ind_age6:ind_age10)/sum(MUemp_age(i,ind_age6:ind_age10))));
    conc_age_11_15(i) = mean(sum(eqbm.conc_age(ind_age11:ind_age15)'.*MUemp_age(i,ind_age11:ind_age15)/sum(MUemp_age(i,ind_age11:ind_age15))));
    conc_age_16_20(i) = mean(sum(eqbm.conc_age(ind_age16:ind_age20)'.*MUemp_age(i,ind_age16:ind_age20)/sum(MUemp_age(i,ind_age16:ind_age20))));
    conc_age_21_25(i) = mean(sum(eqbm.conc_age(ind_age21:ind_age25)'.*MUemp_age(i,ind_age21:ind_age25)/sum(MUemp_age(i,ind_age21:ind_age25))));

    conc_age_above25(i) = mean(sum(eqbm.conc_age(27:end)'.*MUemp_age(i,27:end)/sum(MUemp_age(i,27:end))));

    aggrev(i) = sum(MU_h(i,:).*rev_firms + MU_l(i,:).*rev_firms);
    vw_firms_l(i,:) = MU_l(i,:).*rev_firms/aggrev(i) ;
    vw_firms_h(i,:) = MU_h(i,:).*rev_firms/aggrev(i) ;
    lshare_firms_l = transpose(eqbm.nstar_l)./rev_firms ; % labor share by firm size
    lshare_firms_h = transpose(eqbm.nstar_h)./rev_firms ; % labor share by firm size
    agglshare(i) = sum(lshare_firms_h.*vw_firms_h(i,:))+ sum(lshare_firms_l.*vw_firms_l(i,:)) ;

end

tab = table(year, lf_growth, Ns, ...
    mt, Mt, Xt, E, AFSt, startup_rate, exit_rate, ...
     agglshare, frac250plus_emp, frac10kplus_emp,...
    frac1to4,frac1kplus_emp,avg_fsize_age_6_10,avg_fsize_age_11_15,avg_fsize_age_16_20,...
    avg_fsize_age_21_25,avg_fsize_age_above25,exit_rate_age_6_10,...
    exit_rate_age_11_15,exit_rate_age_16_20,exit_rate_age_21_25,exit_rate_age_above25,...
    conc_age_6_10,conc_age_11_15,conc_age_16_20,conc_age_21_25,conc_age_above25);

results.tab = tab ;
results.MU_h = MU_h;
results.MU_l = MU_l;
results.MUemp_h = MUemp_h;
results.MUemp_l = MUemp_l;
results.MU_age = MU_age ;
results.MUemp_age = MUemp_age ;

f = results ;

