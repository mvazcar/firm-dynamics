% Solve the free-entry condition for the wage w, holding the output price p
% fixed as the numeraire (cf. the original zstar_fun/fzero solve for zstar in
% stationary_alt.m). Only z = p/w is pinned down; the initial guess w = p
% corresponds to z = 1.

function w = solve_wage(params)

fun = @(w) entry_residual(w, params) ;
w = fzero(fun, params.p) ;
