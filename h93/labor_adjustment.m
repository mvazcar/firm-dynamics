% =========================================================================
% Labor adjustment (firing) cost of Hopenhayn & Rogerson (1993):
%     g(n, n_lag) = tau * max(0, n_lag - n).
% A linear, one-sided cost of tau per worker fired -- only shrinking the
% workforce (n < n_lag) is costly; hiring is free. tau is denominated in
% WAGES (HR: "a tax equal to 1 year's wages" is tau = 0.2 for a 5-year
% period), so the dollar cost is w*g. Vectorizes over n and n_lag.
% =========================================================================

function g = labor_adjustment(n, n_lag, tau)
g = tau .* max(0, n_lag - n) ;
