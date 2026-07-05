% =========================================================================
% 2-D value function iteration for Hopenhayn & Rogerson (1993).
% State (s_j, n_lag_i); choice n'_l on the employment grid params.n. Bellman:
%
%   W(s,n_lag) = max_{n'} { return_fn(s,n_lag,n') + beta*max(E_sW(s',n'), -w*tau*n') }
%
% where the nested max is next period's exit decision (stay vs fire-all).
% Returns the value W (nz x na), the employment policy index Npol, the
% survival indicator survX(j,l) (=1 if a firm with shock j and chosen
% employment l stays next period), and the entrant value We(s) (n_lag = 0,
% no cf). Prices only enter through z = p/w, so free entry pins the ratio.
%
% The flow return is split into a state-and-choice part PR(j,l) (revenue,
% wage bill, fixed cost = return_fn with n_lag=0) minus the firing cost
% FC(i,l) = w*labor_adjustment(n(l),n(i),tau), which depends on (n_lag,n').
% =========================================================================

function [W, Npol, survX, We] = vfn_2d(tau, w, params)
nz = params.nz ; nn = params.na ; p = params.p ; theta = params.theta ; beta = params.beta ; n = params.n ;

PR      = p*params.rev0 - w*n' - p*params.cf ;      % nz x na, revenue - wage bill - p*cf
FC      = w*tau*params.fire0 ;                  % na x na, w*labor_adjustment(n(l),n(i),tau)
exitval = -w*tau*n ;                       % na x 1, exit value = -w*g(0,n') = -w*tau*n'

W = zeros(nz,nn) ; Npol = ones(nz,nn) ;
for iter = 1:params.maxit
    EW = params.Pi*W ;                          % nz x na, E_{s'|j}[W(s',l)]
    CV = beta*max(EW, exitval') ;          % nz x na, continuation (stay vs exit)
    A  = PR + CV ;                         % nz x na
    Wn = -inf(nz,nn) ; Np = ones(nz,nn) ;
    for l = 1:nn                           % maximize over the employment choice
        cand = A(:,l) - FC(:,l)' ;         % nz x na: (j,i) = A(j,l) - FC(i,l)
        b = cand > Wn ; Wn(b) = cand(b) ; Np(b) = l ;
    end
    if max(abs(Wn(:)-W(:))) < params.tol, W = Wn ; Npol = Np ; break ; end
    W = Wn ; Npol = Np ;
end

EW    = params.Pi*W ;
survX = EW >= exitval' ;                    % nz x na, stay next period iff EW >= exit
CV    = beta*max(EW, exitval') ;
We    = max(p*params.rev0 - w*n' + CV, [], 2) ;  % nz x 1, entrant value (no cf, n_lag=0)
