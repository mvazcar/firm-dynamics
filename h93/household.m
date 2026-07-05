% =========================================================================
% Representative-household closure for Hopenhayn & Rogerson (1993), with
% indivisible labor (Hansen 1985 / Rogerson 1988):
%     preferences  sum_t beta^t [ ln(c_t) - A * N_t ],   N = fraction employed.
%
% The firm block determines the stationary distribution up to a scale factor
% (the mass of entrants M); the household pins that scale via labor-market
% clearing, which gives HR's employment and welfare numbers. All firm-block
% quantities passed in are PER UNIT of entry mass (M = 1).
%
% Static stationary household problem (1/(1+r) = beta):
%     max ln(c) - A*N   s.t.   p*c = w*N + Pi + R
% FOC gives  c = w/(A*p).  The firing-tax revenue R is rebated lump sum, so it
% cancels in the resource constraint and aggregate consumption equals net
% output:  c = M*c1,  c1 = Y1 - cf*Minc1 - ce*Ment1.  With labor clearing
% N = M*Ld1, the scale is  M = w/(A*p*c1).  A is calibrated at tau = 0 so that
% N = Ntarget (HR use 0.6).
%
% Returns H with the entrant mass M, aggregate employment N, consumption cons,
% labor productivity prod, utility U, and the tau=0-indexed series HR report:
% emp_idx, cons_idx, and uac_idx (utility-adjusted consumption).
% =========================================================================

function H = household(R, P, Ntarget)
p  = P.p ;
c1 = @(r) r.Y - P.cf*r.Minc - P.ce*r.Ment ;      % net output per unit entry mass

% Calibrate A at tau = 0 so that N(0) = Ntarget:
%   N = w*Ld1/(A*p*c1) = Ntarget  ->  A = w0*Ld1(0)/(Ntarget*p*c1(0)).
r0  = R{1} ;
H.A = r0.w*r0.Ld / (Ntarget * p * c1(r0)) ;

nT = numel(R) ;
H.M = zeros(1,nT) ; H.N = H.M ; H.cons = H.M ; H.prod = H.M ; H.U = H.M ;
for it = 1:nT
    r  = R{it} ;  cc = c1(r) ;
    M  = r.w/(H.A*p*cc) ;                          % entry mass (scale)
    H.M(it)    = M ;
    H.N(it)    = M*r.Ld ;                          % aggregate employment
    H.cons(it) = M*cc ;                            % = w/(A*p)
    H.prod(it) = r.Y/r.Ld ;                        % labor productivity (scale-free)
    H.U(it)    = log(H.cons(it)) - H.A*H.N(it) ;
end

% Indices relative to tau = 0 (as reported in HR Table 3)
H.emp_idx  = 100*H.N   ./ H.N(1) ;
H.cons_idx = 100*H.cons./ H.cons(1) ;
H.uac_idx  = 100*exp(H.U - H.U(1)) ;              % utility-adjusted consumption
