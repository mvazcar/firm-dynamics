% =========================================================================
% Maximized per-period profit, in DOLLARS, for the firm problem in
% Hopenhayn, Neira and Singhania (2022), evaluated at the optimal labor
% demand n*(s) from n.m:
%
%   Pi(s) = p*exp(s)*n*(s)^alpha  -  w*n*(s)  -  w*cf
%         = revenue           -  wage bill  -  fixed-cost bill
%
% Both prices appear EXPLICITLY: the output price p multiplies revenue, the
% wage w multiplies the labor bill and the fixed-cost bill (cf is denominated
% in units of labor, so its dollar cost is w*cf).
%
% Equivalently, substituting the first-order condition w*n* = alpha*p*exp(s)*n*^alpha,
%   Pi(s) = (1-alpha)*p*exp(s)*n*(s)^alpha - w*cf.
%
% The value function (vfn) and the entry cost condition are handled in these
% same dollar units; see entry_residual.m (free entry: E[V(s0)] = w*ce).
% =========================================================================

function f = prof_fn(params, s, p, w)
alpha = params.alpha ;
cf    = params.cf ;

nstar = n(params, s, p, w) ;
f = p .* exp(s) .* nstar.^alpha  -  w .* nstar  -  w .* cf ;
