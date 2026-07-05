% =========================================================================
% Solve for the equilibrium wage w that clears free entry, given a hardcoded
% entry cost ce (in goods) and the output price p (the numeraire), for the
% single-type Hopenhayn (1992) model.
%
% Free entry (in dollars) requires  E[V(s0;p,w)] = p*ce  (see entry_residual).
% Only the ratio z = p/w is pinned down, so fixing p and solving for w is a
% normalization choice. The residual E[V(s0)] - p*ce goes from positive (low
% w, high relative price) to negative (high w), so it has a sign change; we
% bracket w from a wide bracket on z = p/w, widen it if necessary, and call
% fzero.
% =========================================================================

function w = solve_wage(params)
p = params.p ;

% Bracket the wage from a wide bracket on the relative price z = p/w:
%   z in [z_lo, z_hi]   <=>   w = p./z in [p/z_hi, p/z_lo].
z_lo = 0.25 ;
z_hi = 4 ;
w_bracket = sort([p/z_hi, p/z_lo]) ;

res_lo = entry_residual(w_bracket(1), params) ;
res_hi = entry_residual(w_bracket(2), params) ;

% Widen the bracket geometrically until it straddles a root.
tries = 0 ;
while sign(res_lo) == sign(res_hi) && tries < 20
    z_lo = z_lo/2 ;
    z_hi = z_hi*2 ;
    w_bracket = sort([p/z_hi, p/z_lo]) ;
    res_lo = entry_residual(w_bracket(1), params) ;
    res_hi = entry_residual(w_bracket(2), params) ;
    tries  = tries + 1 ;
end

if sign(res_lo) == sign(res_hi)
    error('solve_wage:noBracket', ...
        'Could not bracket the free-entry wage; check ce and the price p.') ;
end

w = fzero(@(x) entry_residual(x, params), w_bracket) ;

% Confirm the returned wage actually clears free entry.
res = entry_residual(w, params) ;
if abs(res) > 1e-6 * max(1, abs(p*params.ce))
    error('solve_wage:notCleared', ...
        'Free entry not cleared at w = %.6g (residual = %.3e); check ce.', w, res) ;
end
