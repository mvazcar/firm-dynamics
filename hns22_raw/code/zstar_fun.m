% This code is used in the alternative explanations code, but not in the main code.
% It is used to solve for zstar such that the entry condition holds.

function entryvalue = zstar_fun(zstar,params)
ce = params.ce ;
G = params.G ;
F_h = params.F_h ;
F_l = params.F_l ;
cf_h = params.cf_h;
cf_l = params.cf_l;

omega_h = params.omega_h;

omega_l = 1 - omega_h;

v_h = vfn(F_h, zstar, cf_h, params)  ;

v_l = vfn(F_l, zstar, cf_l, params)  ;

entryvalue = omega_h*sum(v_h.*G) + omega_l*sum(v_l.*G) - ce  ;  
