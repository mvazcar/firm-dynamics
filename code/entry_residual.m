% Free-entry residual as a function of the wage w (price p fixed). Analogous
% to the original zstar_fun.m, but solving for w rather than zstar: it returns
% the expected value of an entrant minus the entry cost. Here firm value is in
% dollars and ce is in units of labor, so the entry cost is w*ce and free
% entry is  omega_h*sum(v_h.*G) + omega_l*sum(v_l.*G) = w*ce.

function res = entry_residual(w, params)

p = params.p ;
G = params.G ;
omega_h = params.omega_h ;
omega_l = 1 - omega_h ;

v_h = vfn(params.F_h, p, w, params.cf_h, params) ;
v_l = vfn(params.F_l, p, w, params.cf_l, params) ;

res = omega_h*sum(v_h.*G) + omega_l*sum(v_l.*G) - w*params.ce ;
