% This function returns the maximized profit function 
% \pi(s,p) = max_n{p*exp(s)*n^alpha - n - cf}

function f = prof_fn(params, s, z)

alpha = params.alpha ;
cf = params.cf ;

prof = (z*exp(s)*alpha^alpha).^(1/(1-alpha))*(1 - alpha) - cf ;
f = prof ;