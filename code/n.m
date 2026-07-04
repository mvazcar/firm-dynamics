% =========================================================================
% Optimal (profit-maximizing) labor demand n*(s) for a firm with
% productivity state s, facing output price p and wage w, in Hopenhayn,
% Neira and Singhania (2022).
%
% FIRM'S STATIC PROBLEM (in dollars):
%       max_n  p*exp(s)*n^alpha  -  w*n  -  w*cf
%
% FIRST-ORDER CONDITION with respect to n:
%       d/dn [ p*exp(s)*n^alpha - w*n - w*cf ] = 0
%       alpha*p*exp(s)*n^(alpha-1) - w = 0
%       n^(1-alpha) = alpha*p*exp(s) / w
% =>    n*(s) = ( alpha*p*exp(s) / w )^(1/(1-alpha)).
%
% Both the output price p and the wage w appear EXPLICITLY: revenue is scaled
% by p, the marginal cost of labor is w. (The fixed cost w*cf is a lump-sum
% overhead and does not enter the interior FOC, so it is absent from n*.)
% =========================================================================

function f = n(params, s, p, w)
alpha = params.alpha ;

f = (alpha .* p .* exp(s) ./ w).^(1/(1-alpha)) ;
