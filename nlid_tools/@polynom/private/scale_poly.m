function newC = scale_poly(C,M);
% shifts and scales a (usually normalized) polynomial.
% syntax:  newC = scale_poly(C,M);
%
%  C represents an N'th order polynomial, and is a N+1 element
%  vector, stored in order of descending powers.  Thus, the polynomial is
%  given by:  C(1)*X^N + ... + C(N)*X + C(N+1).  (see polyval) 
%  M is a 2 element (i.e. first-order polynomial) vector.
%
%  newC returns the coefficients of the polynomial C, evaluated at
%  M(1) X + M(2).

% Copyright 1998-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

C = C(:);
N1 = length(C);
order = N1 - 1;
T = zeros(N1,N1);

% T is a matrix whose columns contain the (first order) polynomial M, raised
% to powers [order ... 1 0 ]
T(N1,N1) = 1;
T(N1-1:N1,N1-1) = M;
for i = 2:order
  temp = conv(T(:,N1-1),T(:,N1-i+1));
  T(:,N1-i) = temp(N1:2*N1-1);
end

newC = T*C;

