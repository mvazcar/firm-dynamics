% =========================================================================
% Value function iteration for an incumbent firm in Hopenhayn, Neira and
% Singhania (2022). Given the output price p, the wage w, the operating cost
% cf and the type-specific transition matrix F, it returns the firm value
% function V in DOLLARS (the same units as prof_fn). The firm exits whenever
% its continuation value would be negative:
%
%       V(s) = max{ 0, prof_fn(s;p,w) + beta * E[V(s')|s] }.
%
% CONVERGENCE NOTE: the original HNS22 code caps the iteration at maxiter=100,
% which does NOT reach tol=1e-8 (the operator is a beta-contraction and needs
% several hundred iterations). We preserve maxiter=100 as the DEFAULT so the
% recode reproduces the published, calibrated benchmark exactly; the cap is
% now exposed through params.maxiter so counterfactuals can request a fully
% converged solve (params.maxiter = 5000, say). See the README.
% =========================================================================

function f = vfn(F, p, w, cf, params)
beta = params.beta ;
svec = params.svec ;
ns   = params.ns ;
tol  = params.tol ;

if isfield(params, 'maxiter')
    maxiter = params.maxiter ;
else
    maxiter = 100 ;            % HNS22 default (see convergence note above)
end

params.cf = cf ;

v  = zeros(ns,1) ;
Tv = v + 1 ;
iter = 0 ;
while max(abs(Tv - v)) > tol && iter <= maxiter
    v = Tv ;
    Ev = transpose(v'*F) ;
    Tv = max(0, prof_fn(params, svec, p, w) + beta*Ev) ;
    iter = iter + 1 ;
end

f = v ;
