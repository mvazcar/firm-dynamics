% =========================================================================
% Stationary equilibrium and aggregates of the 2-D HR1993 firm block, given
% the firing tax tau and the wage w (with p fixed). Builds the policy-induced
% transition over the (s, n_lag) state, solves for the stationary measure
% mu = (I - T)^{-1} (M * entry), and returns the Table-3 firm statistics.
%
% Timing: a firm active at (s_j, n_lag_i) chooses n' = P.n(Npol(j,i)); at the
% start of next period it exits (before the shock) if survX(j,Npol(j,i))=0,
% else it draws s' ~ F(j,.) and is active at (s', n'). Entrants (mass M=1 per
% period) arrive at n_lag = 0 with s ~ v.
% =========================================================================

function r = stationary_2d(tau, w, P)
nz = P.nz ; nn = P.na ; n = P.n ; Pi = P.Pi ; theta = P.theta ; s = P.s ;
[~, Npol, survX, ~] = vfn_2d(tau, w, P) ;

emp = n(Npol) ;                                 % nz x na, employment n'(j,i)
lag = repmat(n', nz, 1) ;                        % nz x na, lagged employment n(i)
survpol = zeros(nz,nn) ;                          % survX(j, Npol(j,i))
for i = 1:nn, for j = 1:nz, survpol(j,i) = survX(j, Npol(j,i)) ; end, end

% ---- Stationary measure over (s, n_lag) --------------------------------
nS = nz*nn ; i0 = 1 ;                             % entrants at n_lag index 1 (n = 0)
II = zeros(nz*nS,1) ; JJ = II ; VV = II ; c = 0 ;
for i = 1:nn
    for j = 1:nz
        if survpol(j,i)
            lp   = Npol(j,i) ;
            rows = (1:nz)' + (lp-1)*nz ;          % destinations (k, lp)
            col  = j + (i-1)*nz ;
            II(c+1:c+nz) = rows ; JJ(c+1:c+nz) = col ; VV(c+1:c+nz) = Pi(j,:)' ;
            c = c + nz ;
        end
    end
end
T   = sparse(II(1:c), JJ(1:c), VV(1:c), nS, nS) ;
psi = zeros(nS,1) ; psi((1:nz)' + (i0-1)*nz) = P.v ;   % entry inflow (M = 1)
mu  = (speye(nS) - T) \ psi ;
MU  = reshape(max(mu,0), nz, nn) ;

% ---- Aggregates (per unit entry mass, M = 1) ---------------------------
Mass = sum(MU(:)) ;  N = sum(sum(MU.*emp)) ;
r.avg_size  = N/Mass ;
r.exit_rate = sum(sum(MU.*(~survpol)))/Mass ;

% Levels needed for the household closure (all per unit of entry mass):
r.Mass = Mass ;                                  % mass of active firms
r.Ld   = N ;                                     % labor demand (employment)
r.Y    = sum(sum(MU.*(s.*(emp.^theta)))) ;       % gross output
r.Ment = sum(P.v) ;                              % mass of entrants (= 1)
r.Minc = Mass - r.Ment ;                         % incumbents (pay cf; entrants do not)

hire = max(0, emp - lag) ; fire = max(0, lag - emp) ;
JC = sum(sum(MU.*hire)) ;                          % entrants counted (n_lag = 0)
JD = sum(sum(MU.*fire)) + sum(sum(MU.*(~survpol).*emp)) ;   % + exiters fire all
r.turnover = (JC+JD)/(2*N) ;

layoff = sum(sum(MU.*(w*tau*fire))) + sum(sum(MU.*(~survpol).*(w*tau*emp))) ;
r.layoff_wb = layoff/(w*N) ;

r.prod = sum(sum(MU.*(s.*(emp.^theta)))) / N ;     % average labor productivity

% ---- Serial correlation and growth variance of log n (survivors) -------
ln = log(max(emp, 1e-12)) ;
Sw=0; Ex=0; Ey=0; Exx=0; Eyy=0; Exy=0 ;
for i = 1:nn
    for j = 1:nz
        if survpol(j,i) && MU(j,i) > 0
            lp = Npol(j,i) ; x = ln(j,i) ;
            wgt = MU(j,i)*Pi(j,:)' ; y = ln(:,lp) ; sw = sum(wgt) ;
            Sw=Sw+sw; Ex=Ex+x*sw; Ey=Ey+wgt'*y;
            Exx=Exx+x^2*sw; Eyy=Eyy+wgt'*(y.^2); Exy=Exy+x*(wgt'*y);
        end
    end
end
Ex=Ex/Sw; Ey=Ey/Sw; Vx=Exx/Sw-Ex^2; Vy=Eyy/Sw-Ey^2; Cxy=Exy/Sw-Ex*Ey;
r.serial_n   = Cxy/sqrt(Vx*Vy) ;
r.var_growth = Vx + Vy - 2*Cxy ;

% ---- Raw objects for plotting (makefigs) -------------------------------
r.MU = MU ; r.Npol = Npol ; r.emp = emp ; r.lag = lag ;
r.survpol = survpol ; r.T = T ; r.n = n ; r.s = s ;
