function Mnew = poly_rescale(Mold,PolyInput)
% transforms polynomial coefficients for different input stats
%
%  syntax:  Mnew = poly_rescape(Mold,PolyInput)
%
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


% first, we extract the details from the polynomial object
C = Mold.polyCoef;
OldType = get(Mold,'polyType');
Range = get(Mold,'polyRange');
xmax = Range(2);
xmin = Range(1);
oldRange = (xmax-xmin);
oldAvg = (xmax+xmin)/2;
oldMean = get(Mold,'polyMean');
oldStd = get(Mold,'polyStd');

newMin = min(PolyInput);
newMax = max(PolyInput);
newRange = newMax - newMin;
newAvg = (newMax + newMin)/2;
newMean = mean(PolyInput);
newStd = std(PolyInput);

% first, we convert the polynomial to a normalized power form
C = flipud(C(:));
N1 = length(C);
order = N1 - 1;

switch lower(OldType(1));
  case 'p'
    % power series -- do nothing
    Cnew = C;   
  case 't'
    % tchebyshev polynomial, relevant stats are  average and range,
    m1 = newRange/oldRange;
    m0 = 2 * (newAvg - oldAvg)/oldRange;
    M = [m1; m0];
    S = scale_matrix(M,order);
    T = tcheb_matrix(order);
    Cnew = T\(S*T*C);
    % compute new values for mean and std.
    newStd = oldStd*m1;
    newMean = (oldMean - oldAvg) *m1 + newAvg;  
  case 'h'
    % hermite polynomial, relevant stats are mean and std
    m1 = newStd/oldStd;
    m0 = (newMean - oldMean)/oldStd;
    M = [m1; m0];
    T = hermite_matrix(order);
    S = scale_matrix(M.order);
    Cnew = T\(S*T*C);
    % compute new values for average and range.
    newRange = oldRange*m1;
    newAvg = (oldAvg - oldMean)*m1 + newMean;
    newMax = newAvg + newRange/2;
    newMin = newAvg - newRange/2;    
  otherwise
    error('OldType not recognized');
end

Cnew = flipud(Cnew);
newRange = [newMin; newMax];

Mnew = Mold;
set(Mnew,'polyMean',newMean,'polyStd',newStd,'polyRange',newRange);
set(Mnew,'polyCoef',Cnew);
return



function T = scale_matrix(M,order);
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
if order==0,
    error ('Order must be > 0');
end

N1 = order+1;
T = zeros(N1,N1);

% T is a matrix whose columns contain the (first order) polynomial M, raised
% to powers [order ... 1 0 ]
T(N1,N1) = 1;
T(N1-1:N1,N1-1) = M;
for i = 2:order
  temp = conv(T(:,N1-1),T(:,N1-i+1));
  T(:,N1-i) = temp(N1:2*N1-1);
end



function T = tcheb_matrix (order);
% computes a matrix that transforms tchebyshev polynomials to power series
% form. 
% syntax:  T = tcheb_matrix (order);
% 
% where order is the order of the polynomial to be transformed, and T is the
% transformation matrix -- then Coeff_power = T Coeff_tchebyshev
% 
% note that coefficients are assumed to be stored in order of decreasing
% exponent (see polyval).

N1 = order+1;
T = zeros(N1,N1);
% set up the zero and first order tchebyshev polynomials. (1 and x)
T(N1,N1) = 1;
T(N1-1,N1-1) = 1;

% use the recurrance relation  T(n+1) = 2x T(n) - T(n-1) 
% to generate the remaining columns.
% David being clever.
for i = 2:order
  T(1:order,N1-i) = 2*T(2:N1,N1-i+1);
  T(:,N1-i) =   T(:,N1-i) - T(:,N1-i+2); 
end

function T = hermite_matrix (order);
% computes a matrix that transforms hermite polynomials to power series
% form. 
% syntax:  T = hermite_matrix (order);
% 
% where order is the order of the polynomial to be transformed, and T is the
% transformation matrix -- then Coeff_power = T Coeff_hermite
% 
% note that coefficients are assumed to be stored in order of decreasing
% exponent (see polyval).

N1 = order+1;
T = zeros(N1,N1);
% set up the zero and first order hermite polynomials. (1 and x)
T(N1,N1) = 1;
T(N1-1,N1-1) = 1;


% use the recurrance relation T(n+1) = xT(n) - (n-1)T(n-1)
% to generate the remaining columns.
for i = 2:order
  T(1:order,N1-i) = T(2:N1,N1-i+1);
  T(:,N1-i) =   T(:,N1-i) - (i-1)*T(:,N1-i+2); 
end












