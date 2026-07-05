% =========================================================================
% Solve for the equilibrium wage w that clears free entry in the 2-D HR1993
% model, given the firing tax tau and the output price p (numeraire).
% Free entry:  sum_s v(s) * W_e(s; p, w) = p*ce.  Because firm value is
% homogeneous of degree 1 in (p,w), this pins only z = p/w; with p fixed we
% solve for w. The entry value is decreasing in w, so we bracket and use fzero.
% =========================================================================

function w = solve_wage_2d(tau, P)
f = @(x) entry_resid(x, tau, P) ;
wlo = P.p/3 ; whi = P.p/0.5 ;
flo = f(wlo) ; fhi = f(whi) ; tries = 0 ;
while sign(flo)==sign(fhi) && tries < 25
    wlo = wlo/1.4 ; whi = whi*1.4 ;
    flo = f(wlo) ; fhi = f(whi) ; tries = tries + 1 ;
end
if sign(flo)==sign(fhi)
    error('solve_wage_2d:noBracket', 'Could not bracket the free-entry wage (tau=%.3g).', tau) ;
end
w = fzero(f, [wlo whi]) ;
end

function r = entry_resid(w, tau, P)
    [~,~,~,We] = vfn_2d(tau, w, P) ;
    r = sum(P.v.*We) - P.p*P.ce ;
end
