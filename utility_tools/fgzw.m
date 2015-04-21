function [c,dcda] = fgzw(t,params)
% compute 2-nd order IRF for k,z,w
% [c dcda] = fkzw(t,a)
% fgzw computes a 2'nd order  impulse response and its partial
% derivatives at each time t.  The derivatives are computed with respect
% to gain (g) , natural frequency (Wn), damping factor (zeta), 
%
% t = vector of times at which the irf is to be evaluated
% param = parameters = [g  zeta Wn]
%             g    = steady state gain
%             zeta = damping parameter
%             Wn   = natural freuqency (rad/s);
%
% Transfer function is:
%
%                g
%          ---------------------
%           
%          s^2     + 2 z     +  1
%          __         __
%          wn^2       wn
%
% nb: partials were computeed using Maple.
% C =  irf computed at each time T
% dcda = partial derivatives w.r.t A at each time T
%
% Caveat Emptor:  the special case of zeta=1 is not handle properly yet
%

% Copyright 1995, Robert E Kearney and James P Trainor
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% See also:  NLLS
% JPT
% 20 Feb 1995 REK develop from fkzw
%
g = params(1);
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

dcda(1,:)= wn*c3.*c1./c4;

dcda(2,:) = -g*wn^2*t.*c1.*c3./c4 + g*z*wn^2*c3.*c2.*t./c4^2 ...
    - g*z*wn*c3.*c1./c4^3;

dcda(3,:) = g*c3.*c1./c4 - g*wn*z*t.*c3.*c1./c4 + g*wn*c3.*c2.*t;


dcda = real(dcda)';

end

%  

c = g* wn*c3.*c1 ./ (c4);


% The partial derivative equations generate "intermediate" complex results
% which combine to give a real results, but due to round off error, the
% derivatives are often of the form: real + 0.000000i.  The identification
% function (nlls) works fine even when the derivatives have this form,
% (I have verified this)  but since I know the result is always real, I will
% due away will the extraneous imaginary part, also, I reshape dcda so that
% it works properly with nlls:
% reshape this too
c = c(:);

% end of file fkzw.m













