% =========================================================================
% Maximized per-period profit, in DOLLARS, for the single-type Hopenhayn
% (1992) firm, evaluated at the optimal labor demand n*(s) from n.m.
%
% COST CONVENTION (this model): the operating cost cf is denominated in
% GOODS (units of the output good), so its dollar cost is p*cf (not w*cf):
%
%   Pi(s) = p*exp(s)*n*(s)^alpha  -  w*n*(s)  -  p*cf
%         = revenue           -  wage bill  -  fixed-cost bill (in goods)
%
% Both prices appear explicitly: p multiplies revenue AND the goods-
% denominated fixed cost; w multiplies the labor bill. Because cf does not
% depend on n, it does not enter the labor first-order condition, so n*(s)
% is unchanged from the standard case (see n.m).
% =========================================================================

function f = prof_fn(params, s, p, w)
alpha = params.alpha ;
cf    = params.cf ;

nstar = n(params, s, p, w) ;
f = p .* exp(s) .* nstar.^alpha  -  w .* nstar  -  p .* cf ;
