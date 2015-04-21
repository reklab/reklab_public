function [c,dcda] = fkzw(t,params)
% compute 2-nd order IRF for k,z,w
% [c dcda] = fkzw(t,a)
% FKZW computes a 2'nd order compliance impulse response and its partial
% derivatives at each time T.  The derivatives are computed with respect
% to natural frequency (Wn), damping factor (zeta), and stiffnes (k).  
%
% T = vector of times at which the irf is to be evaluated
% A = parameters = [k  zeta Wn]
% k    = steady state gain
% zeta = damping parameter
% Wn   = natural freuqency (rad/s);


% Copyright 1995, Robert E Kearney and James P Trainor
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


% C = compliance irf computed at each time T
% DCDA = partial derivatives w.r.t A at each time T
%
% Caveat Emptor:  the special case of zeta=1 is not handle properly yet
%
% See also:  NLLS
% JPT
% 14 Feb 1995 REK Change order of parameters
%                 k is now steady state gain
%
k = params(1);
z = params(2);
wn = params(3);


if z == 1
  error('sorry, this routine can''t yet handle zeta=1')
end

t=[t(:)]';

c1 = sinh(wn*sqrt(-1+z^2)*t);
c2 = cosh(wn*sqrt(-1+z^2)*t);
c3 = exp(-z*wn*t);
c4 = sqrt(-1+z^2);

if nargout == 2

% Note, the derivatives have been generated using Maple.

% > simplify(diff(invlaplace((1/k)/(s^2/(wn^2)+s*2*z/wn+1),s,t),wn));
%                                    2 1/2                             2 1/2
% - exp(- z wn t) (- sinh(wn (- 1 + z )    t) + wn z t sinh(wn (- 1 + z )    t)
% 
%                           2 1/2              2 1/2    /            2 1/2
%      - wn cosh(wn (- 1 + z )    t) t (- 1 + z )   )  /  (k (- 1 + z )   )
%                                                     /

dcda(3,:) = -c3 .* (-c1 + wn*z*c1.*t - wn*c2*c4.*t) / (k*c4);

% > simplify(diff(invlaplace((1/k)/(s^2/(wn^2)+s*2*z/wn+1),s,t),z)); 
%                                                       2 1/2
%            - wn exp(- z wn t) (- wn t sinh(wn (- 1 + z )    t)
%
%                                        2 1/2     2
%                 + wn t sinh(wn (- 1 + z )    t) z
%
%                                      2 1/2                2 1/2
%                 - wn cosh(wn (- 1 + z )    t) t z (- 1 + z )
%
%                                   2 1/2         /            2 3/2
%                 + sinh(wn (- 1 + z )    t) z)  /  (k (- 1 + z )   )
%                                               /

dcda(2,:) = -wn*c3 .* (-wn*t.*c1 + wn*t.*c1*z^2 - wn*c2.*t*z*c4 + c1*z) / ...
            (k*c4^3);

% >  simplify(diff(invlaplace((1/k)/(s^2/(wn^2)+s*2*z/wn+1),s,t),k));
%                                                     2 1/2
%                    wn exp(- z wn t) sinh(wn (- 1 + z )    t)
%                  - -----------------------------------------
%                                  2         2 1/2
%                                 k  (- 1 + z )

dcda(1,:) = -wn*c3.*c1/(k^2*c4);


% The partial derivative equations generate "intermediate" complex results
% which combine to give a real results, but due to round off error, the
% derivatives are often of the form: real + 0.000000i.  The identification
% function (nlls) works fine even when the derivatives have this form,
% (I have verified this)  but since I know the result is always real, I will
% due away will the extraneous imaginary part, also, I reshape dcda so that
% it works properly with nlls:
dcda = real(dcda)';

end

% and the actual impulse response:
% > invlaplace((1/k)/(s^2/(wn^2) + s*2*z/wn + 1),s,t);
%                                                    2 1/2
%                   wn exp(- z wn t) sinh(wn (- 1 + z )    t)
%                   -----------------------------------------
%                                          2 1/2
%                                k (- 1 + z )

c = wn*c3.*c1 / (k*c4);

% reshape this too
c = c(:);

% end of file fkzw.m













