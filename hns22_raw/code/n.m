% ##########################################################################
% Maximized employment function used in Hopenhayn, Neira and Singhania (2022)
% ##########################################################################

function f = n(params, s, z)
alpha = params.alpha ;
f = (z*exp(s)*alpha).^(1/(1-alpha)) ;
