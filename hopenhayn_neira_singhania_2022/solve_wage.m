% Solve the free-entry condition for the wage w, holding the output price p
% fixed as the numeraire. This mirrors the original zstar_fun/fzero solve for
% zstar (stationary_alt.m): only z = p/w is pinned down, so we fix p and use
% fzero on the entry residual. The initial guess w = p corresponds to z = 1.

function w = solve_wage(params)

fun = @(w) entry_residual(w, params) ;
w = fzero(fun, params.p) ;
