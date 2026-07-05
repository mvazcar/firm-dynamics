% =========================================================================
% Free-entry residual as a function of the wage w, with the output price p
% held fixed (the numeraire), for the single-type Hopenhayn (1992) model.
%
% A potential entrant draws its initial productivity from G and enters if the
% expected DOLLAR value of entering covers the DOLLAR entry cost. The entry
% cost ce is denominated in GOODS, so its dollar cost is p*ce, and the free-
% entry (zero expected profit of entry) condition is
%
%   sum( V(s;p,w).*G )  =  p*ce.
%
% entry_residual returns  E[V(s0;p,w)] - p*ce.  The market-clearing wage is
% its root (see solve_wage). Firm value V is homogeneous of degree 1 in
% (p,w) (cf and ce are in goods, so they scale with p), so this condition
% pins down only the ratio z = p/w; with p fixed we solve it for the wage w.
% =========================================================================

function res = entry_residual(w, params)
p = params.p ;
G = params.G ;

V = vfn(params.F, p, w, params.cf, params) ;

entry_value = sum(V.*G) ;              % expected dollar value of entry
res = entry_value - p*params.ce ;      % free entry: E[V(s0)] = p*ce (ce in goods)
