% Solve the free-entry condition for the wage w in the 2-D HR1993 model, with
% the output price p fixed as the numeraire (cf. the original zstar_fun/fzero
% solve for zstar). Free entry:  sum_s v(s) * W_e(s;p,w) = p*ce. Only z = p/w
% is pinned down; the initial guess w = p corresponds to z = 1.

function w = solve_wage_2d(tau, params)

fun = @(w) entry_resid(w, tau, params) ;
w = fzero(fun, params.p) ;
end

function r = entry_resid(w, tau, params)
[~,~,~,We] = vfn_2d(tau, w, params) ;
r = sum(params.v.*We) - params.p*params.ce ;
end
