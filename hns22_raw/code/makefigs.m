% =========================================================================
% This program makes the plots for the firm entry-exit model in Hopenhayn, 
% Neira and Singhania (2022). It loads the calibrated.mat output file, reads 
% data from .csv files stored in the 'data/Calibration' folder, and saves the 
% figures in .eps format in the ../../figures/ folder.
% =========================================================================
clc
clear
close all

set(0,'defaulttextInterpreter','latex') 
saveplots = 1; % 0: don't save plots; 1: save plots for slides (with title); 2: save plots for paper (no title); 3: backup figures

load benchmark

%% Prepare data summary stats
datadir = '../data_summary_stats/' ;

% Read startup rate and other variables from the data so we can compare
startup_data_1940_1962 = csvread(strcat(datadir, 'startup_rate_1940_1960.csv'), 1, 0) ;
startup_data_1980on = csvread(strcat(datadir, 'startup_rate.csv'), 1, 0) ;
exitrate_data_1980on = csvread(strcat(datadir, 'exit_rate.csv'), 1, 0) ;
avgfsize = csvread(strcat(datadir, 'avgfsize.csv'), 1, 0) ;

% Firm-Size and Firm-Emp Distribution 
FSizeEmpDist1978 = csvread(strcat(datadir, 'fsizeempdist.csv'), 1, 0) ;
fsizedist1978 = FSizeEmpDist1978(:,2);
fempdist1978 = FSizeEmpDist1978(:,3);

% Age Distribution: Firm and Employment
AgeDist = csvread(strcat(datadir, 'fagedist.csv'), 1, 0) ;
agefirmdist2014 = AgeDist(:,4);

% New firms to employment ratio
enttoemp_data = csvread(strcat(datadir, 'entranttoemp.csv'), 1, 0) ;
enttoemp_data_years = enttoemp_data(:,1);
enttoemp_data_vals = enttoemp_data(:,2);

shareemp_250plus_data = csvread(strcat(datadir, 'shareemp_250plus.csv'), 1, 0) ;
laborshare_data = csvread(strcat(datadir, 'laborshare.csv'), 1, 0) ;
laborshare_KSZ = csvread(strcat(datadir, 'laborshare_KSZ.csv'), 1, 0) ;

ind1978 = find(bm.year == 1978);
%% Figure 1: Civilian labor force growth rate
% First calculate average labor force growth rate for 10 year intervals
subperiods = transpose(1940:10:2010);

for i=1:length(subperiods)-1
   ind = find(bm.year == subperiods(i)); 
   indprime = find(bm.year == subperiods(i+1)); 
   lfgrowth_agg(i,1) = mean(bm.lf_growth(ind:indprime-1)-1) ;
end

lfgrowth_agg(i+1,1) = mean(bm.lf_growth(indprime:end)-1) ;
lfgrowth_agg = lfgrowth_agg*100 ;

Decade = num2str(subperiods) ;
Decade2 = num2str(subperiods(2:end)) ;


figure(1) % Civilian Labor Force Growth Rate By Decade
bar(subperiods, lfgrowth_agg, 1, 'FaceColor', [0.9290    0.6940    0.1250], 'EdgeColor', [0 0 0], 'LineWidth', 1)
set(gca, 'YGrid', 'on', 'XGrid', 'off')
text(subperiods, lfgrowth_agg, num2str(lfgrowth_agg, '%0.2f'), 'HorizontalAlignment', 'center', 'VerticalAlignment','bottom', 'FontSize', 16); 
ax = gca ;
ax.YLim = [0 3] ;
ax.XTickLabels = cellstr(strcat(Decade, 's')) ;
ax.XAxis.TickLength = [0,0];
ytickformat('%.1f%%')
set(gca,'FontSize',14)
if saveplots == 1
    title('Civilian Labor Force Growth Rate')
    text(1922, -0.3,'Source: BLS Current Population Survey','clipping','off');
    print (gcf,'-depsc2','../../figures/Figure_1')
end
hold off

%% Figure 4: Entry Rate (Data Only)

figure(4) % Entry Rate (Data Only)
hold on
grid on
plot([startup_data_1940_1962(1:end,1); startup_data_1980on(2:end,1)], [startup_data_1940_1962(1:end,2); startup_data_1980on(2:end,2).*100],'-.', 'LineWidth', 4)
plot(startup_data_1940_1962(1:end,1), startup_data_1940_1962(1:end,2),'Color',[0    0.4470    0.7410], 'LineWidth', 4)
plot(startup_data_1980on(2:end,1), startup_data_1980on(2:end,2).*100,'Color',[0    0.4470    0.7410], 'LineWidth', 4)
xlim([1940,2020])
set(gca,'FontSize',16)
ytickformat('%.0f%%')
if saveplots == 1
title('Entry Rate')
        print (gcf,'-depsc2','../../figures/Figure_4')
end
hold off

%% Figure 5: Firm Size Distributions

ind_1978 = find(bm.year == 1978);

Size_Dist_1978Trans = [sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 4)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 4));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 9)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 9))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 4)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 4)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 19)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 19))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 9)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 9)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 49)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 49))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 19)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 19)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 99)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 99))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 49)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 49)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 249)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 249))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 99)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 99)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 499)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 499))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 249)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 249)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 999))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 499)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 499)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 2499)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 2499))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 999)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 4999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 4999))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 2499)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 2499)));
    (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 9999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 9999))) - (sum(results.MU_h(ind_1978,stateqbm.nstar_h <= 4999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l <= 4999)));
    sum(results.MU_h(ind_1978,stateqbm.nstar_h > 9999)) + sum(results.MU_l(ind_1978,stateqbm.nstar_l > 9999))]/(sum(results.MU_l(ind_1978,:))+sum(results.MU_h(ind_1978,:)));

Emp_Dist_1978Trans = [sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 4)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 4));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 9)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 9))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 4)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 4)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 19)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 19))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 9)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 9)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 49)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 49))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 19)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 19)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 99)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 99))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 49)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 49)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 249)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 249))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 99)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 99)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 499)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 499))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 249)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 249)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 999))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 499)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 499)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 2499)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 2499))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 999)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 4999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 4999))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 2499)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 2499)));
    (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 9999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 9999))) - (sum(results.MUemp_h(ind_1978,stateqbm.nstar_h <= 4999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l <= 4999)));
    sum(results.MUemp_h(ind_1978,stateqbm.nstar_h > 9999)) + sum(results.MUemp_l(ind_1978,stateqbm.nstar_l > 9999))]/(sum(results.MUemp_l(ind_1978,:))+sum(results.MUemp_h(ind_1978,:)));


c = categorical({'1 to 4','5 to 9','10 to 19','20 to 49','50 to 99','100 to 249','250 to 499','500 to 999','1000 to 2499','2500 to 4999','5000 to 9999','10,000+'});
c = reordercats(c,{'1 to 4','5 to 9','10 to 19','20 to 49','50 to 99','100 to 249','250 to 499','500 to 999','1000 to 2499','2500 to 4999','5000 to 9999','10,000+'});

figure(5) % 1978 Firm Distributions 
    pos1 = get(gcf,'Position'); 
    set(gcf,'units','points','position',pos1.*[1,1,2,1])

subplot(1, 2, 1) % 1978 Firm Size Distributions
b1 = bar(c,[fsizedist1978 Size_Dist_1978Trans].*100,1);
set(b1(2),'FaceAlpha',0.7);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
legend('Data','Model')
ax = gca ;
ax.XAxis.TickLength = [0,0];
ytickformat('%.0f%%')
title('(a) 1978 Firm Size Distribution')
ylabel('Share of Firms')

subplot(1, 2, 2) % 1978 Firm Employment Distributions
b1 = bar(c,[fempdist1978 Emp_Dist_1978Trans].*100,1);
set(b1(2),'FaceAlpha',0.7);
set(gca, 'YGrid', 'on', 'XGrid', 'off')
legend('Data','Model','Location','northwest')
ax = gca ;
ax.XAxis.TickLength = [0,0];
ytickformat('%.0f%%')
title('(b) 1978 Firm Employment Distribution')
ylabel('Share of Employment')
print(gcf,'-depsc2','../../figures/Figure_5')


%% Figure 6: Entry Rate (Model, Data, Labor Force Growth)

figure(6) % Entry Rate
hold on
grid on
plot([startup_data_1940_1962(1:end,1); startup_data_1980on(2:end,1)], [startup_data_1940_1962(1:end,2); startup_data_1980on(2:end,2).*100],'-.', 'LineWidth', 4)
plot(startup_data_1940_1962(1:end,1), startup_data_1940_1962(1:end,2),'Color',[0    0.4470    0.7410], 'LineWidth', 4)
plot(startup_data_1980on(2:end,1), startup_data_1980on(2:end,2).*100,'Color',[0    0.4470    0.7410], 'LineWidth', 4)
plot(bm.year, bm.startup_rate.*100,'Color',[0.8500    0.3250    0.0980], 'LineWidth', 4)
plot(bm.year, (bm.lf_growth-1).*100,'Color',[0.9290    0.6940    0.1250], 'LineWidth', 4)
xlim([1940,2020])
ylim([-2,18])
htext1 = text(bm.year(40),bm.startup_rate(60).*100-1.7,'Entry Rate, Model');
set(htext1,'Color',[0.8500    0.3250    0.0980],'VerticalAlignment','bottom','FontSize',16) 
htext2 = text(1985,13.2,'Entry Rate, Data');
set(htext2,'Color',[0    0.4470    0.7410],'Rotation',0,'VerticalAlignment','bottom','FontSize',16)
htext3 = text(1978.5,2.8,'Labor force growth rate, Data');
set(htext3,'Color',[0.9290    0.6940    0.1250],'Rotation',0,'VerticalAlignment','bottom','FontSize',16)
set(gca,'FontSize',14)
ytickformat('%.0f%%')
title('Entry Rate')
hold off
print(gcf,'-depsc2','../../figures/Figure_6')


%% Figure 7: Firm Age Distributions

bm_agedist =results.MU_age./sum(results.MU_age,2);
norise_agedist =results_norise.MU_age./sum(results_norise.MU_age,2);
notrans_agedist =results_notrans.MU_age./sum(results_notrans.MU_age,2);

bm_ageempdist =results.MUemp_age./sum(results.MUemp_age,2);
norise_ageempdist =results_norise.MUemp_age./sum(results_norise.MUemp_age,2);
notrans_ageempdist =results_notrans.MUemp_age./sum(results_notrans.MUemp_age,2);

ind1977=find(bm.year==1977);


figure(7) 
    pos1 = get(gcf,'Position'); 
    set(gcf,'units','points','position',pos1.*[1,1,2,1])

subplot(1, 2, 1) % Firm age distribution 1978 - CDF
hold on
grid on
plot(stateqbm.age,cumsum(bm_agedist(ind1977+1,:)).*100,'Color',[0.8500    0.3250    0.0980], 'LineWidth', 4);
plot(stateqbm.age,cumsum(norise_agedist(ind1977+1,:)).*100, '-.','Color',[0.9290    0.6940    0.1250],'LineWidth', 4);
xlim([0 30])
legend('Transition','Steady State','Location','SouthEast')
set(gca,'FontSize',14)
ytickformat('%.0f%%')
title('(a) 1978 CDF')



subplot(1, 2, 2) % Firm age distribution 2014 - CDF
hold on
grid on
h1 = plot([0;5;10;15;20;25], (1-agefirmdist2014)*100, '.', 'Color', [0    0.4470    0.7410], 'MarkerSize', 30);
plot(stateqbm.age,cumsum(bm_agedist(end,:)).*100, 'LineWidth', 4);
plot(stateqbm.age,cumsum(notrans_agedist(end,:)).*100, '-.','LineWidth', 4);
uistack(h1,'top')
xlim([0 30])
ylim([0 100])
legend('Transition','Steady State','Data','Location','SouthEast')
set(gca,'FontSize',14)
ytickformat('%.0f%%')
title('(b) 2014 CDF')
print(gcf,'-depsc2','../../figures/Figure_7')


%% Figure 8: Aggregate Exit, Average Firm Size and Concentration Time Series

figure(8) % Exit, Size, and Concentration 
    pos1 = get(gcf,'Position'); 
    set(gcf,'units','points','position',pos1.*[1,1,2,1])

subplot(1, 3, 1) % Aggregate Exit Rate since 1978
hold on
grid on
plot(exitrate_data_1980on(1:end,1), exitrate_data_1980on(1:end,2).*100,'Color',[0    0.4470    0.7410], 'LineWidth', 4)
plot(1978:2014,exitrate_data_1980on(1,2).*100 - 0.054*(0:36),'-.','Color',[0.9290    0.6940    0.1250],'LineWidth', 4)
plot(bm.year, bm.exit_rate.*100,'Color',[0.8500    0.3250    0.0980  0.7], 'LineWidth', 4)
xlim([1978,2020])
ylim([7,11.5])
set(gca,'FontSize',14)
ytickformat('%.0f%%')
pos1 = get(gcf,'Position'); 
set(gcf,'Position', pos1 - [pos1(3),0,0,0]) 
title('(b) Exit Rate')
legend('Data','Data, Aging Only','Model','Location','SouthWest')
hold off


subplot(1, 3, 2) % Average firm size since 1978
hold on
grid on
plot(avgfsize(2:end-1, 1), avgfsize(2:end-1, 2), 'LineWidth', 4)
plot(1978:2014,avgfsize(2, 2)+0.239*(0:36),'-.','Color',[0.9290    0.6940    0.1250],'LineWidth', 4)
plot(bm.year, bm.AFSt,'Color',[0.8500    0.3250    0.0980  0.7],'LineWidth', 4)
xlim([1978,2020])
pos3 = get(gcf,'Position');  
set(gcf,'Position', pos3 + [pos3(3),0,0,0]) 
set(gca,'FontSize',14)
title('(a) Average firm size')
legend('Data','Data, Aging Only','Model','Location','SouthEast')
hold off



subplot(1, 3, 3) % Concentration since 1978
hold on
grid on
plot(shareemp_250plus_data(:,1),shareemp_250plus_data(:,2).*100, 'LineWidth', 4);
plot(1978:2014,shareemp_250plus_data(1,2).*100 + 0.228*(0:36),'-.','Color',[0.9290    0.6940    0.1250],'LineWidth', 4)
plot(bm.year,bm.frac250plus_emp.*100,'Color',[0.8500    0.3250    0.0980 0.7], 'LineWidth', 4);
xlim([1978,2020])
set(gca,'FontSize',14)
ytickformat('%.0f%%')
pos12 = get(gcf,'Position');   
set(gcf,'Position', pos12 + [pos12(3)/2,0,0,0]) 
title('(c) Concentration')
legend('Data','Data, Aging Only','Model','Location','SouthEast')
hold off

print(gcf,'-depsc2','../../figures/Figure_8')


%% Figure 9: Corporate Labor Share

figure(9) % Labor Share, Koh Santaeulalia-Llopis Zheng 
hold on
grid on
h1 = plot(laborshare_data(:,1),(laborshare_data(:,2)-laborshare_data(6,2)).*100, '-.', 'LineWidth', 4);
plot(laborshare_KSZ(:,1),(laborshare_KSZ(:,2)-laborshare_KSZ(33,2)).*100,'Color', [0.9290    0.6940    0.1250], 'LineWidth', 4);
plot(bm.year,(bm.agglshare-bm.agglshare(42)).*100, 'Color',[0.8500    0.3250    0.0980],'LineWidth', 4)
xlim([1940,2020])
legend('Karabarbounis and Neiman (2014)','Koh, Santaeulalia-Llopis and Zheng (2020)','Model','Location','SouthWest')
uistack(h1,'top')
set(gca,'FontSize',14)
ytickformat('%.0f%%')
pos13 = get(gcf,'Position');   
set(gcf,'Position', pos13 + [pos13(3)/2,0,0,0]) 
title('Corporate Labor Share')
hold off
print(gcf,'-depsc2','../../figures/Figure_9')


%% Figure 11: Growth in the number of firms

enttoemp = bm.mt./bm.E;

figure(11) % Ratio of new firms to total employment
hold on
grid on
plot(enttoemp_data_years, enttoemp_data_vals*100, 'LineWidth', 4)
plot(enttoemp_data_years, enttoemp(ind1978:end)*100, 'LineWidth', 4)
xlim([1978,2020])
pos3 = get(gcf,'Position');  
set(gcf,'Position', pos3 + [pos3(3),0,0,0]) 
set(gca,'FontSize',14)
ytickformat('%.2f%%')
title('Ratio of New Firms to Total Employment')
hold off
print(gcf,'-depsc2','../../figures/Figure_11')
