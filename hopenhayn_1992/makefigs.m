% =========================================================================
% Figures for the single-type Hopenhayn (1992) model (hopenhayn_1992), on the shared
% HR1993 calibration. Saves PNGs to hopenhayn_1992/figures/.
%   Fig 1  Firm and employment size distributions (model vs HR data).
%   Fig 2  Value function V(s), labor demand n*(s), and the exit threshold.
%   Fig 3  Cohort survival and exit hazard by age.
% =========================================================================
clc ; clear ; close all ;
outdir = fullfile(fileparts(mfilename('fullpath')), 'figures') ;
if ~exist(outdir,'dir'), mkdir(outdir) ; end

% ---- Build the model (same calibration as main.m) ----------------------
p.beta=0.8; p.alpha=0.64; p.ns=100; p.tol=1e-8; p.maxiter=1000; p.mstar=1;
p.ce=14.9087; p.cf=15.1537; p.p=1;
rho=0.93; sigma_eps=(1-p.alpha)*sqrt(0.53); mu=0.8707; nmax=5000; gridR=5.5; vfrac=0.74;
logSmax=(1-p.alpha)*log(nmax)-log(p.alpha);
svec=linspace(logSmax-gridR,logSmax,p.ns)';
p.svec=svec; p.F=transpose(tauchen_fixed(svec,mu,rho,sigma_eps));
nlow=floor(vfrac*p.ns); G=zeros(p.ns,1); G(1:nlow)=1/nlow; p.G=G;

w=solve_wage(p); p.w=w;
eq=stationary(p);
V=vfn(p.F, p.p, p.w, p.cf, p);          % value function (hopenhayn_1992 timing)
nstar=eq.nstar; mustar=eq.mustar; sstar=eq.sstar_ind;
s=exp(svec);

col=[0 0.30 0.65];

% ---- Fig 1: size distributions -----------------------------------------
edges=[0 20 100 500 inf]; lbl={'1-19','20-99','100-499','500+'};
fs=zeros(1,4); es=zeros(1,4);
for b=1:4, in=nstar>=edges(b)&nstar<edges(b+1); fs(b)=sum(mustar(in)); es(b)=sum(mustar(in).*nstar(in)); end
fs=fs/sum(fs); es=es/sum(es);
f1=figure('Position',[100 100 820 380],'Color','w');
subplot(1,2,1); bar([fs; .52 .37 .10 .01]'); set(gca,'XTickLabel',lbl,'FontSize',11);
title('Firm size distribution'); ylabel('share of firms'); legend({'model','HR data'});
subplot(1,2,2); bar([es; .06 .24 .37 .33]'); set(gca,'XTickLabel',lbl,'FontSize',11);
title('Employment size distribution'); ylabel('share of employment'); legend({'model','HR data'});
savefig_png(f1, fullfile(outdir,'fig1_size_distribution.png'));

% ---- Fig 2: value function and labor demand ----------------------------
f2=figure('Position',[100 100 820 400],'Color','w');
subplot(1,2,1); hold on;
plot(svec, V, '-','Color',col,'LineWidth',2);
yl=ylim; plot(svec(sstar)*[1 1], yl, 'k--');
grid on; set(gca,'FontSize',11); xlabel('log productivity, log s'); ylabel('value V(s)');
title('Value function'); legend({'V(s)','exit threshold s^*'},'Location','northwest');
subplot(1,2,2);
semilogy(svec, nstar, '-','Color',col,'LineWidth',2); hold on;
yl=ylim; plot(svec(sstar)*[1 1], yl, 'k--');
grid on; set(gca,'FontSize',11); xlabel('log productivity, log s'); ylabel('employment n^*(s)');
title('Optimal labor demand'); legend({'n^*(s)','exit threshold s^*'},'Location','northwest');
savefig_png(f2, fullfile(outdir,'fig2_value_policy.png'));

% ---- Fig 3: cohort survival and hazard ---------------------------------
Gnew=G; Gnew(1:sstar-1)=0; Fnew=p.F; Fnew(:,1:sstar-1)=0;
c=Gnew; S=zeros(1,16); S(1)=sum(c);
for t=1:15, c=Fnew*c; S(t+1)=sum(c); end
surv=S/S(1); haz=(S(1:end-1)-S(2:end))./S(1:end-1);
f3=figure('Position',[100 100 820 400],'Color','w');
subplot(1,2,1);
plot(0:15, surv,'-o','Color',col,'LineWidth',2,'MarkerFaceColor',col); grid on;
set(gca,'FontSize',11); xlabel('cohort age (periods)'); ylabel('survival share');
title('Cohort survival');
subplot(1,2,2);
plot(1:15, haz,'-o','Color',col,'LineWidth',2,'MarkerFaceColor',col); grid on;
set(gca,'FontSize',11); xlabel('cohort age (periods)'); ylabel('exit hazard');
title('Cohort exit hazard');
savefig_png(f3, fullfile(outdir,'fig3_cohort_survival.png'));

fprintf('Saved 3 figures to %s\n', outdir);

function savefig_png(f, path)
    try, exportgraphics(f, path, 'Resolution', 130);
    catch, print(f, path, '-dpng', '-r130'); end
end
