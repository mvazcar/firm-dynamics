% =========================================================================
% Hopenhayn & Rogerson (1993) BENCHMARK (tau = 0), firm block only.
% Calibrates the model to the paper's targets and reports Table 2, using the
% paper's structural parameters and the correct HR exit timing.
%
% HR timing (paper p.7): the exit decision is made at the start of a period
% BEFORE the new shock, so the value recursion puts the max on the
% CONTINUATION:
%     W(s) = pi_gross(s) - p*cf + beta*max( E_{s'|s} W(s'), 0 )      (tau = 0)
% with pi_gross(s) = max_n { p*s*n^theta - w*n } and s = exp(log s) in LEVELS.
% A firm with current shock s exits next period iff E_s W(s') < 0 (s <= s*).
% Entrants draw s ~ v, pay NO cf (footnote 5): W_e(s) = W(s) + p*cf, and free
% entry sets ce = E_v[W_e]/p.
%
% Normalization p = w = 1 (HR benchmark). Grid: fixed, uniform in log s, top
% state set so n*(s_max) = 5000 (HR). Calibration:
%   - vfrac (entrant dist = uniform over bottom vfrac of the grid) pinned so
%     that the average entrant size = 7.5;
%   - (mu_s, cf) calibrated so that average firm size = 61.7 and the 5-year
%     exit rate = 0.37.
% =========================================================================
clc ; clear ; close all ;

% ---- Structural parameters (paper) -------------------------------------
P.beta    = 0.8 ;                     % 5-year discount factor
P.theta   = 0.64 ;                    % production curvature f = s*n^theta
P.rho     = 0.93 ;                    % persistence of log s
P.sigma   = (1-P.theta)*sqrt(0.53) ;  % innovation std of log s (= 0.2621)
P.ns      = 100 ;                     % productivity grid points
P.nmax    = 5000 ;                    % top state implies ~5000 employees (HR)
P.gridR   = 5.5 ;                     % log-grid range below the top
P.p       = 1 ; P.w = 1 ;             % benchmark normalization
P.tol     = 1e-10 ; P.maxit = 5000 ;

% ---- Targets (Table 1 / 2) ---------------------------------------------
TGT.size = 61.7 ;   % mean employment
TGT.exit = 0.37 ;   % 5-year exit rate
TGT.new  = 7.5  ;   % average size of a new firm

% ---- Fixed grid (independent of the calibration) -----------------------
theta = P.theta ;
logSmax = (1-theta)*log(P.nmax) - log(theta) ;   % n*(s_max) = nmax at p=w=1
G.logs  = linspace(logSmax - P.gridR, logSmax, P.ns)' ;
G.s     = exp(G.logs) ;
G.nstar = (theta*P.p.*G.s/P.w).^(1/(1-theta)) ;  % static labor demand

% ---- Pin vfrac so the average entrant size hits the target -------------
% avg_new = mean(nstar(1:nlow)) is monotone in nlow and grid-only.
nlow = find(cumsum(G.nstar)'./(1:P.ns) >= TGT.new, 1) ;
if isempty(nlow), nlow = P.ns ; end
vfrac = nlow / P.ns ;
G.v = zeros(P.ns,1) ; G.v(1:nlow) = 1/nlow ;

% ---- Calibrate (mu_s, cf) to (average size, exit rate) -----------------
[mu_s, cf] = calibrate(TGT, G, P) ;
st = solve_bench(mu_s, cf, G, P) ;

% ---- Report vs Table 2 -------------------------------------------------
fprintf('\n===== HR1993 BENCHMARK (tau = 0) : calibration =====\n') ;
fprintf('  sigma_eps = %.4f   rho = %.2f   beta = %.2f   theta = %.2f\n', P.sigma, P.rho, P.beta, P.theta) ;
fprintf('  calibrated mu_s (mean log s) = %.4f  ->  a = %.4f\n', mu_s, (1-P.rho)*mu_s) ;
fprintf('  calibrated cf               = %.4f\n', cf) ;
fprintf('  entrant dist: uniform over bottom %d/%d states (vfrac = %.2f)\n', nlow, P.ns, vfrac) ;
fprintf('  implied entry cost ce       = %.4f\n', st.ce) ;
fprintf('  calibration hit: size %.2f (tgt %.1f), exit %.3f (tgt %.2f), new %.2f (tgt %.1f)\n', ...
        st.avg_size, TGT.size, st.exit_rate, TGT.exit, st.avg_new, TGT.new) ;
fprintf('  grid: n*(s) in [%.3f, %.0f],  s* at index %d/%d\n', G.nstar(1), G.nstar(end), st.sstar_ind, P.ns) ;

fprintf('\n  %-32s %10s %10s\n', 'Statistic', 'Model', 'HR Table2') ;
fprintf('  %-32s %10.2f %10.2f\n', 'Average firm size',        st.avg_size,   61.2) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Co-worker mean',           st.coworker,   747) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Exit rate of firms',       st.exit_rate,  0.39) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Job turnover rate',        st.turnover,   0.30) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Frac hiring by new firms', st.frac_hire_new, 0.15) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Average size of new firm', st.avg_new,    7.5) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Var of growth (survivors)',st.var_growth, 0.55) ;
fprintf('  %-32s %10.2f %10.2f\n', 'Serial corr log n (surv)', st.serial_n,   0.92) ;
fprintf('  %-32s %10.2f %10.2f\n', '  (uncond. autocorr log s)', st.uncond_ac, 0.93) ;

fprintf('\n  Size distribution (firms):   1-19  20-99  100-499  500+\n') ;
fprintf('    Model : %6.2f %6.2f %8.2f %6.2f\n', st.firmshare) ;
fprintf('    HR    : %6.2f %6.2f %8.2f %6.2f\n', [0.52 0.37 0.10 0.01]) ;
fprintf('  Size distribution (employment):\n') ;
fprintf('    Model : %6.2f %6.2f %8.2f %6.2f\n', st.empshare) ;
fprintf('    HR    : %6.2f %6.2f %8.2f %6.2f\n', [0.06 0.24 0.37 0.33]) ;

fprintf('\n  Cohort exit hazard (ages 1,2,5,10):\n') ;
fprintf('    Model : %6.2f %6.2f %6.2f %6.2f\n', st.hazard) ;
fprintf('    HR    : %6.2f %6.2f %6.2f %6.2f\n', [0.75 0.32 0.15 0.10]) ;

save hr_benchmark

% =========================================================================
% Local functions
% =========================================================================
function [mu_s, cf] = calibrate(TGT, G, P)
    cf_lo = 1 ; cf_hi = 80 ;
    for oit = 1:60
        cf = 0.5*(cf_lo + cf_hi) ;
        mu_s = mu_for_size(cf, TGT.size, G, P) ;
        st = solve_bench(mu_s, cf, G, P) ;
        if st.exit_rate > TGT.exit, cf_hi = cf ; else, cf_lo = cf ; end
        if abs(st.exit_rate - TGT.exit) < 1e-4, break ; end
    end
    mu_s = mu_for_size(cf, TGT.size, G, P) ;
end

function mu_s = mu_for_size(cf, target_size, G, P)
    lo = -1 ; hi = 5 ;
    for it = 1:80
        mu_s = 0.5*(lo + hi) ;
        st = solve_bench(mu_s, cf, G, P) ;
        if st.avg_size > target_size, hi = mu_s ; else, lo = mu_s ; end
        if abs(st.avg_size - target_size) < 1e-3, break ; end
    end
end

function st = solve_bench(mu_s, cf, G, P)
    theta = P.theta ; beta = P.beta ; p = P.p ; w = P.w ; ns = P.ns ;
    s = G.s ; nstar = G.nstar ; v = G.v ;

    Pi = tauchen_fixed(G.logs, mu_s, P.rho, P.sigma) ;   % Pr(s'=k | s=j)

    % Unconditional autocorrelation of the discretized log-s chain (diagnostic)
    pis = ones(1,ns)/ns ;
    for it = 1:5000
        pn = pis*Pi ; if max(abs(pn-pis))<1e-12, pis=pn; break ; end ; pis = pn ;
    end
    Els = pis*G.logs ; Vls = pis*(G.logs.^2) - Els^2 ;
    st_uncond_ac = ((pis.*G.logs')*Pi*G.logs - Els^2)/Vls ;

    % Operating profit (gross of cf) and HR-timing value function
    pigross = p*s.*nstar.^theta - w*nstar ;
    W = zeros(ns,1) ;
    for it = 1:P.maxit
        Wn = pigross - p*cf + beta*max(Pi*W, 0) ;
        if max(abs(Wn - W)) < P.tol, W = Wn ; break ; end
        W = Wn ;
    end
    EW   = Pi*W ;
    surv = EW >= 0 ;
    sstar_ind = find(surv, 1) ;

    % Stationary distribution of active firms by current shock (M = 1)
    psi  = (eye(ns) - Pi'*diag(surv)) \ v ;
    psi  = max(psi, 0) ;
    Mass = sum(psi) ;  N = sum(psi.*nstar) ;

    % --- Statistics ---
    st.sstar_ind = sstar_ind ;
    st.uncond_ac = st_uncond_ac ;
    st.ce        = sum(v.*(W + p*cf))/p ;
    st.avg_size  = N / Mass ;
    st.coworker  = sum(psi.*nstar.^2) / N ;
    st.exit_rate = sum(psi(~surv)) / Mass ;
    st.avg_new   = sum(v.*nstar) / sum(v) ;

    exp_up = zeros(ns,1) ; exp_dn = zeros(ns,1) ;
    for j = 1:ns
        exp_up(j) = Pi(j,:) * max(0, nstar - nstar(j)) ;
        exp_dn(j) = Pi(j,:) * max(0, nstar(j) - nstar) ;
    end
    hire_entrants = sum(v.*nstar) ;
    JC = sum(surv.*psi.*exp_up) + hire_entrants ;
    JD = sum(surv.*psi.*exp_dn) + sum((~surv).*psi.*nstar) ;
    st.turnover      = (JC + JD) / (2*N) ;
    st.frac_hire_new = hire_entrants / JC ;

    edges = [0 20 100 500 inf] ; fs = zeros(1,4) ; es = zeros(1,4) ;
    for b = 1:4
        in = nstar >= edges(b) & nstar < edges(b+1) ;
        fs(b) = sum(psi(in)) ; es(b) = sum(psi(in).*nstar(in)) ;
    end
    st.firmshare = fs / sum(fs) ;
    st.empshare  = es / sum(es) ;

    ln = log(nstar) ;
    jm = (surv.*psi) .* Pi ; jm = jm / sum(jm(:)) ;
    ex = sum(sum(jm .* ln')) ; ey = sum(jm,2)' * ln ;
    dln = ln' - ln ;
    st.var_growth = sum(sum(jm .* (dln.^2))) - (sum(sum(jm.*dln)))^2 ;
    vx = sum(jm,1) * (ln.^2) - ex^2 ;
    vy = sum(jm,2)' * (ln.^2) - ey^2 ;
    cxy = sum(sum(jm .* (ln * ln'))) - ex*ey ;
    st.serial_n = cxy / sqrt(vx*vy) ;

    c = v ; S = zeros(1,11) ; S(1) = sum(c) ;
    for t = 1:10
        c = Pi' * (surv.*c) ; S(t+1) = sum(c) ;
    end
    haz = (S(1:end-1) - S(2:end)) ./ S(1:end-1) ;
    st.hazard = haz([1 2 5 10]) ;
end
