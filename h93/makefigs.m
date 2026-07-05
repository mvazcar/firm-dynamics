% =========================================================================
% Figures for the Hopenhayn & Rogerson (1993) replication (h93). The paper
% itself reports results in tables; these plots visualize the same objects.
% Saves PNGs to h93/figures/.
%   Fig 1  Employment decision-rule bands [n_l(s), n_u(s)] vs productivity,
%          for tau = 0, 0.1, 0.2  (the visual form of HR Table 4).
%   Fig 2  Table-3 policy effects vs tau (model vs HR).
%   Fig 3  Firm and employment size distributions (model vs HR data).
%   Fig 4  Cohort exit-hazard by age (model vs HR).
%   Fig 5  Stationary measure over (s, n_lag) at tau = 0.2 (with the bands).
% =========================================================================
clc ; clear ; close all ;
outdir = fullfile(fileparts(mfilename('fullpath')), 'figures') ;
if ~exist(outdir,'dir'), mkdir(outdir) ; end

% ---- Calibration + grids (moderate resolution for figures) -------------
P.beta=0.8; P.theta=0.64; P.rho=0.93; P.sigma=(1-0.64)*sqrt(0.53);
P.mu=0.8707; P.cf=15.1537; P.ce=14.9087; P.p=1; P.vfrac=0.74;
P.nz=60; P.na=200; P.nzR=5.5; P.nmax=5000; P.tol=1e-7; P.maxit=2000;
logSmax=(1-P.theta)*log(P.nmax)-log(P.theta);
P.logs=linspace(logSmax-P.nzR,logSmax,P.nz)'; P.s=exp(P.logs);
P.Pi=tauchen_fixed(P.logs,P.mu,P.rho,P.sigma);
P.n=[0, logspace(log10(1e-3),log10(1.2*P.nmax),P.na-1)]';
P.rev0=P.s*(P.n'.^P.theta); P.fire0=max(0,P.n-P.n');
nlow=floor(P.vfrac*P.nz); P.v=zeros(P.nz,1); P.v(1:nlow)=1/nlow;

taus=[0 0.1 0.2]; R=cell(1,3);
for it=1:3
    w=solve_wage_2d(taus(it),P);
    R{it}=stationary_2d(taus(it),w,P); R{it}.w=w; R{it}.tau=taus(it);
end
H=household(R,P,0.6);
col=[0 0.30 0.65; 0.85 0.33 0.10; 0.30 0.65 0.20];   % tau colors

% ---- Fig 1: decision-rule bands ----------------------------------------
NL=nan(P.nz,3); NU=nan(P.nz,3); NPT=nan(P.nz,3);
for it=1:3
    Np=R{it}.Npol;
    for j=1:P.nz
        inact=find(Np(j,:)==(1:P.na));            % n' = n_lag  (inaction)
        NPT(j,it)=P.n(Np(j,1));                   % policy at n_lag=0 (hire target)
        if ~isempty(inact), NL(j,it)=P.n(inact(1)); NU(j,it)=P.n(inact(end)); end
    end
end
f1=figure('Position',[100 100 1080 460],'Color','w');
% Left: the bands (zoomed to the mass of firms, log scale)
subplot(1,2,1); hold on;
plot(P.logs, NPT(:,1), '-', 'Color',col(1,:), 'LineWidth',2) ;             % tau=0 line
for it=2:3
    plot(P.logs, NL(:,it), '-',  'Color',col(it,:), 'LineWidth',1.8) ;
    plot(P.logs, NU(:,it), '--', 'Color',col(it,:), 'LineWidth',1.8) ;
end
set(gca,'YScale','log','FontSize',12); grid on;
xlabel('log productivity, log s'); ylabel('employment (log scale)');
title('Decision-rule bands [n_l(s), n_u(s)]');
legend({'\tau=0 (n^*)','\tau=0.1 lower','\tau=0.1 upper','\tau=0.2 lower','\tau=0.2 upper'}, ...
       'Location','northwest');
xlim([1 3.5]); ylim([5 3e3]);
% Right: fractional band width  (n_u-n_l)/midpoint  -- HR: ~1/3 at tau=0.2
subplot(1,2,2); hold on;
for it=2:3
    wdt=(NU(:,it)-NL(:,it))./((NU(:,it)+NL(:,it))/2);
    plot(P.logs, wdt, '-','Color',col(it,:),'LineWidth',2);
end
grid on; set(gca,'FontSize',12); xlabel('log productivity, log s');
ylabel('band width / midpoint'); title('Inaction-band width (distortion)');
legend({'\tau=0.1','\tau=0.2'},'Location','northeast'); xlim([1 3.5]);
sgtitle('Employment decision rules with firing costs','FontSize',13);
savefig_png(f1, fullfile(outdir,'fig1_decision_bands.png')) ;

% ---- Fig 2: Table-3 effects vs tau -------------------------------------
f2=figure('Position',[100 100 820 620],'Color','w') ;
hr.prod=[100 99.2 97.9]; hr.emp=[100 98.3 97.5]; hr.uac=[100 98.7 97.2]; hr.turn=[.30 .26 .22];
mo.prod=100*H.prod/H.prod(1); mo.emp=H.emp_idx; mo.uac=H.uac_idx;
mo.turn=[R{1}.turnover R{2}.turnover R{3}.turnover];
panels={'Average productivity','prod'; 'Total employment','emp'; ...
        'Utility-adjusted consumption','uac'; 'Job turnover rate','turn'};
for k=1:4
    subplot(2,2,k); hold on;
    plot(taus, mo.(panels{k,2}), '-o','Color',col(2,:),'LineWidth',2,'MarkerFaceColor',col(2,:));
    plot(taus, hr.(panels{k,2}), '--s','Color',[.4 .4 .4],'LineWidth',1.5);
    grid on; set(gca,'FontSize',11); xlabel('firing tax \tau'); xticks(taus);
    title(panels{k,1}); if k==1, legend({'model','HR (1993)'},'Location','southwest'); end
end
sgtitle('HR1993 Table 3: effect of the firing tax','FontSize',13);
savefig_png(f2, fullfile(outdir,'fig2_table3_effects.png')) ;

% ---- Fig 3: size distributions (tau=0) ---------------------------------
edges=[0 20 100 500 inf]; lbl={'1-19','20-99','100-499','500+'};
MU=R{1}.MU; emp=R{1}.emp; fs=zeros(1,4); es=zeros(1,4);
for b=1:4, in=emp>=edges(b)&emp<edges(b+1); fs(b)=sum(MU(in)); es(b)=sum(MU(in).*emp(in)); end
fs=fs/sum(fs); es=es/sum(es);
f3=figure('Position',[100 100 820 380],'Color','w');
subplot(1,2,1); bar([fs; .52 .37 .10 .01]'); set(gca,'XTickLabel',lbl,'FontSize',11);
title('Firm size distribution'); ylabel('share of firms'); legend({'model','HR data'});
subplot(1,2,2); bar([es; .06 .24 .37 .33]'); set(gca,'XTickLabel',lbl,'FontSize',11);
title('Employment size distribution'); ylabel('share of employment'); legend({'model','HR data'});
savefig_png(f3, fullfile(outdir,'fig3_size_distribution.png')) ;

% ---- Fig 4: cohort exit hazard (tau=0) ---------------------------------
i0=1; psi=zeros(P.nz*P.na,1); psi((1:P.nz)'+(i0-1)*P.nz)=P.v;
c=psi; S=zeros(1,16); S(1)=sum(c);
for t=1:15, c=R{1}.T*c; S(t+1)=sum(c); end
haz=(S(1:end-1)-S(2:end))./S(1:end-1);
f4=figure('Position',[100 100 720 460],'Color','w'); hold on;
plot(1:15, haz, '-o','Color',col(1,:),'LineWidth',2,'MarkerFaceColor',col(1,:));
plot([1 2 5 10],[.75 .32 .15 .10],'s','Color',[.85 .1 .1],'MarkerSize',10,'LineWidth',2,'MarkerFaceColor',[.85 .1 .1]);
grid on; set(gca,'FontSize',12); xlabel('cohort age (periods)'); ylabel('exit hazard');
title('Cohort exit hazard by age (\tau=0)'); legend({'model','HR data (ages 1,2,5,10)'});
savefig_png(f4, fullfile(outdir,'fig4_cohort_hazard.png')) ;

% ---- Fig 5: stationary measure over (s, n_lag) at tau=0.2 --------------
MU2=R{3}.MU; dens=MU2; dens(dens<1e-12)=NaN; dens=log10(dens);
f5=figure('Position',[100 100 760 540],'Color','w');
imagesc(P.logs, log10(max(P.n,1e-3)), dens', 'AlphaData',~isnan(dens')); axis xy;
set(gca,'Color',[.95 .95 .95],'FontSize',12); colorbar;
xlabel('log productivity, log s'); ylabel('log lagged employment, log n_{-1}');
title('Stationary measure log_{10}\mu(s, n_{-1}) at \tau=0.2'); ylim([-1 4]);
savefig_png(f5, fullfile(outdir,'fig5_stationary_measure.png')) ;

fprintf('Saved 5 figures to %s\n', outdir) ;

function savefig_png(f, path)
    try, exportgraphics(f, path, 'Resolution', 130) ;
    catch, print(f, path, '-dpng', '-r130') ; end
end
