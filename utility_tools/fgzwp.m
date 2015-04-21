function [h,dcda] = fgzwp(t,params)
% compute 3rd order IRF for k,z,w,p
% [h dcda] = fgzwp(t,a)
% FGZWP computes a 3rd order impulse response and its partial
% derivatives at each time T.  The derivatives are computed with respect
% to natural frequency (Wn), damping factor (zeta), steady state gain (g)
% and real pole (P).
%
% T = vector of times at which the irf is to be evaluated
% A = parameters = [k  zeta Wn P]
% g    = steady state gain
% zeta = damping parameter
% Wn   = natural freuqency (rad/s);
% P = real pole
% C = compliance irf computed at each time T
% DCDA = partial derivatives w.r.t A at each time T
%
% Caveat Emptor:  the special case of zeta=1 is not handle properly yet
%
% See also:  NLLS
% Mehdi
% 28 March 1995 
% 26 Aug 98 REK fixed problem with initialization of dcda
%			   for the case of one output variable. 
%
g = params(1);
z = params(2);
wn = params(3);
p = params(4);

if z == 1
  error('sorry, this routine can''t yet handle zeta=1')
end

t=[t(:)]';


beta=p/(z*wn);
alpha1=exp(-z*wn*t);
alpha3=cos(sqrt(1-z^2)*wn*t);
alpha4=sin(sqrt(1-z^2)*wn*t);
alpha5=exp(-p*t);
alpha6=sqrt(1-z^2);


% since beta*z^2*(beta-2)+1=z^2*(beta-1)^2+(1-z^2)>0;
% Coefficient of exp(-pt) is always negative.
%---------------------------------------------------------
%h(t) = - g wn p (z alpha1 alpha4 wn - alpha1 alpha4 p
%             + alpha1 alpha6 alpha3 wn - alpha5 wn alpha6)

%                /    2                2        2 1/2
%               /  ((p  - 2 p z wn + wn ) (1 - z )   )
%              /

h = -g*wn*p*(z*alpha1.*alpha4*wn-alpha1.*alpha4*p+alpha1*alpha6.*alpha3*wn...
-alpha5*wn*alpha6)/(p^2-2*p*z*wn+wn^2)/sqrt(1-z^2);

%--------------------------------------------------------



if nargout == 2


%---------------------------------------------------------------------
%
%          dh(t)
%          ----- = - wn p (z alpha1 alpha4 wn - alpha1 alpha4 p
%            dg

%               + alpha1 alpha6 alpha3 wn - alpha5 wn alpha6)
%
%                 /    2                2        2 1/2
%                /  ((p  - 2 p z wn + wn ) (1 - z )   )
%               /


dcda(1,:) = -wn*p*(z*alpha1.*alpha4*wn-alpha1.*alpha4*p+alpha1*alpha6.*alpha3*...
wn-alpha5*wn*alpha6)/(p^2-2*p*z*wn+wn^2)/sqrt(1-z^2);
%--------------------------------------------------------------------------

%dh(t)          2                  2  3      3                   2
%----- = (- 2 wn  t alpha1 alpha4 p  z  + 2 z  alpha1 alpha4 p wn
%  dz
%
%         2                   2      2   2                         2
%    - 4 z  alpha1 alpha4 wn p  + 2 z  wn  t alpha1 alpha6 alpha3 p
%
%        3                    2                              2  2
%    + wn  t alpha1 alpha4 p z  + 2 alpha1 alpha6 alpha3 p wn  z
%
%       2                     3                2         2
%    + z  wn t alpha1 alpha4 p  - 2 p alpha5 wn  alpha6 z
%
%                                   3                       2
%    - z wn t alpha1 alpha6 alpha3 p  + alpha1 alpha4 p z wn
%
%            3                                2                  2
%    - 3 z wn  t alpha1 alpha6 alpha3 p + 2 wn  t alpha1 alpha4 p  z
%
%                     3       4
%    + alpha1 alpha4 p  z + wn  t alpha1 alpha6 alpha3
%
%                   2                          3                     2
%    + 2 p alpha5 wn  alpha6 - alpha1 alpha4 wn  + alpha1 alpha4 wn p
%
%        3                                                2
%    - wn  t alpha1 alpha4 p - 2 alpha1 alpha6 alpha3 p wn
%
%                          3     2                         2           /
%    - wn t alpha1 alpha4 p  + wn  t alpha1 alpha6 alpha3 p ) p wn g  /
%                                                                    /
%
%          2 3/2   4      3           2   2      2  2   2           3     4
%   ((1 - z )    (p  - 4 p  z wn + 2 p  wn  + 4 p  z  wn  - 4 p z wn  + wn ))


dcda(2,:) =(-2*wn^2*t.*alpha1.*alpha4*p^2*z^3+2*z^3*alpha1.*alpha4*p*wn^2-4*z ...
^2*alpha1.*alpha4*wn*p^2+2*z^2*wn^2*t.*alpha1*alpha6.*alpha3*p^2+wn^3*t.* ...
alpha1.*alpha4*p*z^2+2*alpha1*alpha6.*alpha3*p*wn^2*z^2+z^2*wn*t.*alpha1.* ...
alpha4*p^3-2*p*alpha5*wn^2*alpha6*z^2-z*wn*t.*alpha1*alpha6.*alpha3*p^3+ ...
alpha1.*alpha4*p*z*wn^2-3*z*wn^3*t.*alpha1*alpha6.*alpha3*p+2*wn^2*t.*alpha1.* ...
alpha4*p^2*z+alpha1.*alpha4*p^3*z+wn^4*t.*alpha1*alpha6.*alpha3+2*p*alpha5*wn ...
^2*alpha6-alpha1.*alpha4*wn^3+alpha1.*alpha4*wn*p^2-wn^3*t.*alpha1.*alpha4*p- ...
2*alpha1*alpha6.*alpha3*p*wn^2-wn*t.*alpha1.*alpha4*p^3+wn^2*t.*alpha1*alpha6.* ...
alpha3*p^2)*p*wn*g/(1-z^2)^(3/2)/(p^4-4*p^3*z*wn+2*p^2*wn^2+4*p^2*z^2*wn^2 ...
-4*p*z*wn^3+wn^4);

%--------------------------------------------------------------------

%  dh(t)                3                                3
%  ----- = g p (alpha1 p  alpha4 - z alpha1 alpha4 t wn p
%   dwn

%                         3                       2  2
% - 3 z alpha1 alpha4 t wn  p + alpha1 alpha4 t wn  p

%                     4      2                                      2
% + alpha1 alpha4 t wn  + 2 p  alpha5 wn alpha6 - alpha1 alpha4 p wn
%
%                              3                            3
% + alpha1 alpha6 alpha3 t wn p  + alpha1 alpha6 alpha3 t wn  p
%
%             2                                 2
% - 2 alpha1 p  alpha6 alpha3 wn - 2 p alpha5 wn  alpha6 z
%
%                2                  2      2   2                  2
% - 2 z alpha1 wn  t alpha6 alpha3 p  + 2 z  wn  t alpha1 alpha4 p
%
%                                2      2                   2
% + 2 z alpha1 alpha6 alpha3 p wn  + 2 z  alpha1 alpha4 p wn
%
%                         2    /
% - 2 z alpha1 alpha4 wn p )  /  (
%                            /
%
%      2 1/2   4      3           2   2      2  2   2           3     4
%(1 - z )    (p  - 4 p  z wn + 2 p  wn  + 4 p  z  wn  - 4 p z wn  + wn )
%)


dcda(3,:) = g*p*(alpha1*p^3.*alpha4-z*alpha1.*alpha4.*t*wn*p^3-3*z*alpha1 ...
.*alpha4.*t*wn^3*p+alpha1.*alpha4.*t*wn^2*p^2+alpha1.*alpha4.*t*wn^4+2*p^2* ...
alpha5*wn*alpha6-alpha1.*alpha4*p*wn^2+alpha1*alpha6.*alpha3.*t*wn*p^3+ ...
alpha1*alpha6.*alpha3.*t*wn^3*p-2*alpha1*p^2*alpha6.*alpha3*wn-2*p*alpha5* ...
wn^2*alpha6*z-2*z*alpha1*wn^2.*t*alpha6.*alpha3*p^2+2*z^2*wn^2*t.*alpha1 ...
.*alpha4*p^2+2*z*alpha1*alpha6.*alpha3*p*wn^2+2*z^2*alpha1.*alpha4*p*wn^2- ...
2*z*alpha1.*alpha4*wn*p^2)/sqrt(1-z^2)/(p^4-4*p^3*z*wn+2*p^2*wn^2+4*p^2 ...
*z^2*wn^2-4*p*z*wn^3+wn^4);


%---------------------------------------------------------------------

%dh(t)         2             2
%----- = - g wn  (- alpha5 wn  alpha6 - 2 alpha1 alpha4 p wn
%  dp
%
%                     2                          2    3
% + z alpha1 alpha4 wn  + alpha1 alpha6 alpha3 wn  + p  t alpha5 alpha6
%
%      2                                     2
% - 2 p  t alpha5 wn alpha6 z + p t alpha5 wn  alpha6
%
%                    2    2                         2
% + z alpha1 alpha4 p  + p  alpha5 alpha6 - alpha1 p  alpha6 alpha3)
%
%   /
%  /  (
% /
%
%      2 1/2   4      3           2   2      2  2   2           3     4
%(1 - z )    (p  - 4 p  z wn + 2 p  wn  + 4 p  z  wn  - 4 p z wn  + wn )
%)

%dcda(4,:) = -g*wn^2*(-alpha5*wn^2*alpha6-2*alpha1.*alpha4*p*wn+z*alpha1.*...
%alpha4*wn^2+alpha1*alpha6.*alpha3*wn^2+p^3*t.*alpha5*alpha6-2*p^2*t.*alpha5*...
%wn*alpha6*z+p*t.*alpha5*wn^2*alpha6+z*alpha1.*alpha4*p^2+p^2*alpha5*alpha6-...
%alpha1*p^2*alpha6.*alpha3)/(1-z^2)^(1/2)/(p^4-4*p^3*z*wn+2*p^2*wn^2+4*p^2*z...
%^2*wn^2-4*p*z*wn^3+wn^4);

dcda(4,:)=g*wn^2*(-p^3*t.*alpha5*alpha6-alpha5*alpha6*p^2-alpha1*alpha6 ...
.*alpha3*wn^2+2*alpha1*p.*alpha4*wn-z*alpha1.*alpha4*wn^2-p*t.*alpha5*wn^2* ...
alpha6+2*p^2*t.*alpha5*wn*alpha6*z+alpha5*wn^2*alpha6-z*alpha1*p^2.* ... 
alpha4+alpha1*p^2*alpha6.*alpha3)/sqrt(1-z^2)/(p^4-4*p^3*z*wn+2*p^2*wn ...
^2+4*p^2*z^2*wn^2-4*p*z*wn^3+wn^4);
dcda =real(dcda)';
end
h = h(:);












