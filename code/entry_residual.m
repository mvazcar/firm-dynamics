% =========================================================================
% Free-entry residual as a function of the wage w, with the output price p
% held fixed (the numeraire). A potential entrant draws its initial
% productivity from G and enters if the expected DOLLAR value of entering
% covers the DOLLAR entry cost. The entry cost is ce units of labor, i.e.
% w*ce dollars, so the free-entry (zero expected profit of entry) condition is
%
%   omega_h * sum( V_h(s;p,w).*G )  +  omega_l * sum( V_l(s;p,w).*G )  =  w*ce.
%
% entry_residual returns  E[V(s0;p,w)] - w*ce.  The market-clearing wage is
% its root (see solve_wage). Because the firm value V is homogeneous of
% degree 1 in (p,w) and ce is in labor units, this condition depends on p and
% w only through the ratio z = p/w; with p fixed we solve it for the wage w.
% =========================================================================

function res = entry_residual(w, params)
p       = params.p ;
G       = params.G ;
omega_h = params.omega_h ;
omega_l = 1 - omega_h ;

V_h = vfn(params.F_h, p, w, params.cf_h, params) ;
V_l = vfn(params.F_l, p, w, params.cf_l, params) ;

entry_value = omega_h*sum(V_h.*G) + omega_l*sum(V_l.*G) ;   % expected dollar value of entry
res = entry_value - w*params.ce ;                           % free entry: E[V(s0)] = w*ce
