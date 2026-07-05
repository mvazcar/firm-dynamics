% Free-entry residual as a function of the wage w (price p fixed), for the
% single-type Hopenhayn (1992) model. Analogous to the original zstar_fun.m,
% solving for w rather than zstar: the expected value of an entrant minus the
% entry cost. ce is in goods, so the entry cost is p*ce and free entry is
%   sum( V(s;p,w).*G ) = p*ce.

function res = entry_residual(w, params)

p = params.p ;
G = params.G ;

V = vfn(params.F, p, w, params.cf, params) ;

res = sum(V.*G) - p*params.ce ;
