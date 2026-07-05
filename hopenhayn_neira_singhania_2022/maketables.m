% =========================================================================
% This program makes the tables for the firm entry-exit model in Hopenhayn, 
% Neira and Singhania (2022). 
% =========================================================================
clc
clear
close all

% Load results for Tables 3, 4 and 5
load benchmark 

% Load results for Table 6: Alternative Experiments
load alteqbm_bench 
load alteqbm_lpop 
load alteqbm_lrow 
load alteqbm_hentry
load alteqbm_hrts

format bank % Display tables to two decimal point precision

%% Parameters and Match

Parameter_assigned = {'Discount factor'; 
        'Curvature';
        'Labor Force Growth Rate (SS)'};

Parameter_estimated = {'Entry Cost'; 
        'Mean of entrant dist'; 
        'Std. dev. of entrant dist';
        'AR(1) Persistence' ;
        'AR(1) std. dev.';
        'Long-run mean low';
        'Long-run mean high';
        'Operating cost low-type';
        'Operating cost high-type';
        'Share of high mean startups'} ;

Value_ass = [beta; alpha; g];

moment = {'Entry Rate 78-83'; 
        'Avg. Entrant Size 78-83'; 
        'Conc. of entrants, 78-83';
        '5-year growth rate, 78-83' ;
        '5-year exit rate, 78-83';
        'Emp Share 10,000+ 78-83';
        'Firm Share 1 to 4 78-83';
        'Emp Share 1000+ 78-83';
        'Agg. Concentration 78-83'} ;

temp1 = cumprod(1 + stateqbm.cond_growth_age(1:5));
Five_year_growth = temp1(end)-1 ;

temp2 = struct2cell(targets) ;
moment_data = transpose(horzcat(temp2{:})) ;

ind1978 = find(bm.year == 1978);
ind1983 = find(bm.year == 1983);

moment_model = [mean(bm.startup_rate(ind1978:ind1983))*100;
    stateqbm.avg_stsize;
    stateqbm.conc_age(1)*100;
    Five_year_growth*100;
    (1-stateqbm.St(6))*100;
    mean(bm.frac10kplus_emp(ind1978:ind1983))*100;
    mean(bm.frac1to4(ind1978:ind1983))*100;
    mean(bm.frac1kplus_emp(ind1978:ind1983))*100;
    mean(bm.frac250plus_emp(ind1978:ind1983))*100];

match = table(moment, moment_data, moment_model) ;

param_ass = table(Parameter_assigned, Value_ass) ;
param_estim = table(Parameter_estimated, Value) ;
disp('--------------------------------------------')
disp('************* Table 3 **********************')
disp('--------------------------------------------')
disp(param_ass)
disp(param_estim)
disp(match)

ind1987 = find(bm.year == 1987);
ind1992 = find(bm.year == 1992);
ind1997 = find(bm.year == 1997);
ind2002 = find(bm.year == 2002);
ind2014 = find(bm.year == 2014);

disp('--------------------------------------------')
disp('************* Table 4 **********************')
disp('--------------------------------------------')
disp('Age Profiles')
disp('--------------------------------------------')
Name = {'0'; '1'; '2'; '3';'4';'5';'6 to 10';'11 to 15';'16 to 20';'21 to 25';'Above 25'};
Exit_Rate = [0 ;
    stateqbm.exit_rate_age(1);
    stateqbm.exit_rate_age(2);
    stateqbm.exit_rate_age(3);
    stateqbm.exit_rate_age(4);
    stateqbm.exit_rate_age(5);
    mean(bm.exit_rate_age_6_10(ind1987:ind2014));
    mean(bm.exit_rate_age_11_15(ind1992:ind2014));
    mean(bm.exit_rate_age_16_20(ind1997:ind2014));
    mean(bm.exit_rate_age_21_25(ind2002:ind2014));
    mean(bm.exit_rate_age_above25(ind2002:ind2014))]*100;

Avg_FSize = [stateqbm.avg_fsize_age(1);
    stateqbm.avg_fsize_age(2);
    stateqbm.avg_fsize_age(3);
    stateqbm.avg_fsize_age(4);
    stateqbm.avg_fsize_age(5);
    stateqbm.avg_fsize_age(6);
    mean(bm.avg_fsize_age_6_10(ind1987:ind2014));
    mean(bm.avg_fsize_age_11_15(ind1992:ind2014));
    mean(bm.avg_fsize_age_16_20(ind1997:ind2014));
    mean(bm.avg_fsize_age_21_25(ind2002:ind2014));
    mean(bm.avg_fsize_age_above25(ind2002:ind2014))];

Concentration =[stateqbm.conc_age(1);
    stateqbm.conc_age(2);
    stateqbm.conc_age(3);
    stateqbm.conc_age(4);
    stateqbm.conc_age(5);
    stateqbm.conc_age(6);
    mean(bm.conc_age_6_10(ind1987:ind2014));
    mean(bm.conc_age_11_15(ind1992:ind2014));
    mean(bm.conc_age_16_20(ind1997:ind2014));
    mean(bm.conc_age_21_25(ind2002:ind2014));
    mean(bm.conc_age_above25(ind2002:ind2014))]*100;

T = table(Exit_Rate,Avg_FSize,Concentration,'RowNames',Name);
disp(T)

%% Entry Rate Decomposition from 1978 to 2014
tab1977ind = find(bm.year == 1977) ;
tab1978ind = find(bm.year == 1978) ;
tab2013ind = find(bm.year == 2013) ;
tab2014ind = find(bm.year == 2014) ;
Res_1978 = (bm.startup_rate(tab1978ind)-(bm.lf_growth(tab1978ind)-1)-bm.exit_rate(tab1978ind)+(bm.AFSt(tab1978ind)/bm.AFSt(tab1977ind)-1))*100;
Res_2014 = (bm.startup_rate(tab2014ind)-(bm.lf_growth(tab2014ind)-1)-bm.exit_rate(tab2014ind)+(bm.AFSt(tab2014ind)/bm.AFSt(tab2013ind)-1))*100;

nr1977ind = find(norise.year == 1977) ;
nr1978ind = find(norise.year == 1978) ;
nt2013ind = find(notrans.year == 2013) ;
nt2014ind = find(notrans.year == 2014) ;
LR_Res_1978 = (norise.startup_rate(nr1978ind)-(norise.lf_growth(nr1978ind)-1)-norise.exit_rate(nr1978ind)+(norise.AFSt(nr1978ind)/norise.AFSt(nr1977ind)-1))*100;
LR_Res_2014 = (notrans.startup_rate(nt2014ind)-(notrans.lf_growth(nt2014ind)-1)-notrans.exit_rate(nt2014ind)+(notrans.AFSt(nt2014ind)/notrans.AFSt(nt2013ind)-1))*100;

Res_1978_SS = (norise.startup_rate(nr1978ind)-(norise.lf_growth(nr1978ind)-1)-norise.exit_rate(nr1978ind)+(norise.AFSt(nr1978ind)/norise.AFSt(nr1977ind)-1))*100;
Res_2014_SS = (notrans.startup_rate(nt2014ind)-(notrans.lf_growth(nt2014ind)-1)-notrans.exit_rate(nt2014ind)+(notrans.AFSt(nt2014ind)/notrans.AFSt(nt2013ind)-1))*100;

disp('--------------------------------------------')
disp('************* Table 5 **********************')
disp('--------------------------------------------')
disp('Decomposition: Steady State vs. Transitional Effects')
disp('--------------------------------------------')
Name = {'LF Growth'; 'Exit Rate'; 'AFS Growth'; 'Residual';'Entry Rate'};
Total_Change = [((bm.lf_growth(tab2014ind)-1)-(bm.lf_growth(tab1978ind)-1))*100;
    (bm.exit_rate(tab2014ind)-bm.exit_rate(tab1978ind))*100;
    ((bm.AFSt(tab2014ind)/bm.AFSt(tab2013ind)-1)-(bm.AFSt(tab1978ind)/bm.AFSt(tab1977ind)-1))*100;
    Res_2014-Res_1978;
    (bm.startup_rate(tab2014ind)-bm.startup_rate(tab1978ind))*100] ;
Long_Run_Effect = [((notrans.lf_growth(nt2014ind)-1)-(norise.lf_growth(nr1978ind)-1))*100;
    (notrans.exit_rate(nt2014ind)-norise.exit_rate(nr1978ind))*100;
    ((notrans.AFSt(nt2014ind)/notrans.AFSt(nt2013ind)-1)-(norise.AFSt(nr1978ind)/norise.AFSt(nr1977ind)-1))*100;
    LR_Res_2014-LR_Res_1978;
    (notrans.startup_rate(nt2014ind)-norise.startup_rate(nr1978ind))*100] ;
Trans_Effect_1978 = [((norise.lf_growth(nr1978ind)-1)-(bm.lf_growth(tab1978ind)-1))*100;
    (norise.exit_rate(nr1978ind)-bm.exit_rate(tab1978ind))*100;
    ((norise.AFSt(nr1978ind)/norise.AFSt(nr1977ind)-1)-(bm.AFSt(tab1978ind)/bm.AFSt(tab1977ind)-1))*100;
    Res_1978_SS-Res_1978;
    (norise.startup_rate(nr1978ind)-bm.startup_rate(tab1978ind))*100] ;
Trans_Effect_2014 = [((bm.lf_growth(tab2014ind)-1)-(notrans.lf_growth(nt2014ind)-1))*100;
    (bm.exit_rate(tab2014ind)-notrans.exit_rate(nt2014ind))*100;
    ((bm.AFSt(tab2014ind)/bm.AFSt(tab2013ind)-1)-(notrans.AFSt(nt2014ind)/notrans.AFSt(nt2013ind)-1))*100;
    Res_2014-Res_2014_SS;
    (bm.startup_rate(tab2014ind)-notrans.startup_rate(nt2014ind))*100] ;

T = table(Total_Change,Long_Run_Effect,Trans_Effect_1978,Trans_Effect_2014,'RowNames',Name);
disp(T)

format short % Return to regular formatting

delta_exit_agecond_lpop = ((sum(alteqbm_lpop.exit_rate_age'.*alteqbm_bench.age_dist(1:end-1)/sum(alteqbm_bench.age_dist(1:end-1))))-alteqbm_bench.exit_rate)*100;
delta_afs_agecond_lpop = sum(alteqbm_lpop.avg_fsize_age.*alteqbm_bench.age_dist) - alteqbm_bench.avg_fsize ;
age_empdist_lpop = (alteqbm_lpop.avg_fsize_age.*alteqbm_bench.age_dist)./sum(alteqbm_lpop.avg_fsize_age.*alteqbm_bench.age_dist);
delta_conc_agecond_lpop = (sum(alteqbm_lpop.conc_age.*age_empdist_lpop) - alteqbm_bench.conc)*100;

delta_exit_agecond_hrts = ((sum(alteqbm_hrts.exit_rate_age'.*alteqbm_bench.age_dist(1:end-1)/sum(alteqbm_bench.age_dist(1:end-1))))-alteqbm_bench.exit_rate )*100;
delta_afs_agecond_hrts =  sum(alteqbm_hrts.avg_fsize_age.*alteqbm_bench.age_dist) - alteqbm_bench.avg_fsize; 
age_empdist_hrts = (alteqbm_hrts.avg_fsize_age.*alteqbm_bench.age_dist)./sum(alteqbm_hrts.avg_fsize_age.*alteqbm_bench.age_dist);
delta_conc_agecond_hrts = (sum(alteqbm_hrts.conc_age.*age_empdist_hrts) - alteqbm_bench.conc)*100;

delta_exit_agecond_hentry = ((sum(alteqbm_hentry.exit_rate_age'.*alteqbm_bench.age_dist(1:end-1)/sum(alteqbm_bench.age_dist(1:end-1))))-alteqbm_bench.exit_rate)*100;
delta_afs_agecond_hentry =  sum(alteqbm_hentry.avg_fsize_age.*alteqbm_bench.age_dist)- alteqbm_bench.avg_fsize;
age_empdist_hentry = (alteqbm_hentry.avg_fsize_age.*alteqbm_bench.age_dist)./sum(alteqbm_hentry.avg_fsize_age.*alteqbm_bench.age_dist);
delta_conc_agecond_hentry = (sum(alteqbm_hentry.conc_age.*age_empdist_hentry) - alteqbm_bench.conc)*100;

delta_exit_agecond_lrow = ((sum(alteqbm_lrow.exit_rate_age'.*alteqbm_bench.age_dist(1:end-1)/sum(alteqbm_bench.age_dist(1:end-1))))-alteqbm_bench.exit_rate )*100;
delta_afs_agecond_lrow =  sum(alteqbm_lrow.avg_fsize_age.*alteqbm_bench.age_dist) - alteqbm_bench.avg_fsize; 
age_empdist_lrow = (alteqbm_lrow.avg_fsize_age.*alteqbm_bench.age_dist)./sum(alteqbm_lrow.avg_fsize_age.*alteqbm_bench.age_dist);
delta_conc_agecond_lrow = (sum(alteqbm_lrow.conc_age.*age_empdist_lrow) - alteqbm_bench.conc)*100;

disp('--------------------------------------------')
disp('************* Table 6 **********************')
disp('--------------------------------------------')
disp('Alternative Explanations')
disp('--------------------------------------------')
Name = {'Data (1978-2014)';'Our model with transition'; 'Our model no transition'; 'Increase in returns to scale'; 'Higher entry cost'; 'Decrease in mean reversion'};
Entry_Rate = [-0.056; bm.startup_rate(tab2014ind)-bm.startup_rate(tab1978ind); alteqbm_lpop.startup_rate-alteqbm_bench.startup_rate; alteqbm_hrts.startup_rate-alteqbm_bench.startup_rate; alteqbm_hentry.startup_rate-alteqbm_bench.startup_rate; alteqbm_lrow.startup_rate-alteqbm_bench.startup_rate].*100 ;
Exit_Rate = [-0.016; bm.exit_rate(tab2014ind)-bm.exit_rate(tab1978ind); alteqbm_lpop.exit_rate-alteqbm_bench.exit_rate; alteqbm_hrts.exit_rate-alteqbm_bench.exit_rate; alteqbm_hentry.exit_rate-alteqbm_bench.exit_rate; alteqbm_lrow.exit_rate-alteqbm_bench.exit_rate].*100 ;
Avg_Size = [3.5; bm.AFSt(tab2014ind)-bm.AFSt(tab1978ind); alteqbm_lpop.avg_fsize-alteqbm_bench.avg_fsize; alteqbm_hrts.avg_fsize-alteqbm_bench.avg_fsize; alteqbm_hentry.avg_fsize-alteqbm_bench.avg_fsize; alteqbm_lrow.avg_fsize-alteqbm_bench.avg_fsize] ;
Conc = [0.057; bm.frac250plus_emp(tab2014ind)-bm.frac250plus_emp(tab1978ind); alteqbm_lpop.conc-alteqbm_bench.conc; alteqbm_hrts.conc-alteqbm_bench.conc; alteqbm_hentry.conc-alteqbm_bench.conc; alteqbm_lrow.conc-alteqbm_bench.conc].*100 ;
AgeCond_Exit = [0.4; 0; delta_exit_agecond_lpop; delta_exit_agecond_hrts; delta_exit_agecond_hentry; delta_exit_agecond_lrow] ;
AgeCond_AFS = [-5.2; 0; delta_afs_agecond_lpop; delta_afs_agecond_hrts; delta_afs_agecond_hentry; delta_afs_agecond_lrow] ;
AgeCond_Conc = [-2.5; 0; delta_conc_agecond_lpop; delta_conc_agecond_hrts; delta_conc_agecond_hentry; delta_conc_agecond_lrow] ;
T = table(Entry_Rate,Exit_Rate,Avg_Size,Conc,AgeCond_Exit,AgeCond_AFS,AgeCond_Conc,'RowNames',Name);
T.Variables =  round(T.Variables,1);
disp(T)
