% This function takes model parameters as inputs, does the value function 
% iteration and returns the value function

function f = vfn(F, zstar, cf, params)
beta = params.beta ;
svec = params.svec ;
ns = params.ns ;
tol = params.tol ;

params.cf = cf;

v = zeros(ns,1) ;
Tv = v + 1 ;
maxiter = 100 ;
iter = 0 ;
while max(abs(Tv - v)) > tol && iter <= maxiter
    v = Tv ;
    Ev = transpose(v'*F) ;
    Tv = max(0, prof_fn(params, svec, zstar) + beta*Ev) ;
    iter = iter + 1 ;
end

f = v ;
