function Cnew = poly_convert(C,InputStats,OldType,NewType)
% converts bewteen power, tchebyshev and hermite polynomials
%
%  syntax:  Cnew = poly_convert(C,InputStats,OldType,NewType)
%
% C is a vector of polynomial coefficients stored in order of increasing
% exponent (UNLIKE polyval)
%  
% InputStats = [Xmax, Xmin, Xmean, Xstd];
% OldType and NewType can be "power", "tchebyshev" or "hermite"

% Copyright 1998-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

% first, we convert the polynomial to a normalized power form

C = flipud(C(:));
N1 = length(C);
order = N1 - 1;
xmax = InputStats(1);
xmin = InputStats(2);
xmean = InputStats(3);
xstd = InputStats(4);
xrng = (xmax-xmin);
xavg = (xmax+xmin)/2;


if strcmp(lower(OldType),'power')
  % power series -- scale polynomial to [-1 1] range
  scale = [xrng/2; xavg];
  coeff = scale_poly(C,scale);
elseif strcmp(lower(OldType),'tcheb')
  T = tcheb_matrix(order);          % Tchebyshev transformation matrix
  coeff = T*C;
elseif strcmp(lower(OldType),'hermite')
  H = hermite_matrix(order);          % Hermite transformation matrix
  Ctemp = H*C;                        % power, assuming Hermite scaling
  gain = xrng/(2*xstd);               % Do the max - min scaling 
  offset = (xavg - xmean)/xstd;       % and undo the mean - std scaling.
  scale = [gain; offset];
  coeff = scale_poly(Ctemp,scale);
else
  error('OldType not recognized');
end
  

% now, convert into the destination type.



if strcmp(lower(NewType),'power')
    % power series -- scale polynomial from [-1 1] to Range
    scale = [2/xrng; -2*xavg/xrng];
    Cnew = scale_poly(coeff,scale);
elseif strcmp(lower(NewType),'tcheb')
    T = tcheb_matrix(order);           % Tchebyshev transformation matrix
    Cnew = T\coeff;
elseif strcmp(lower(NewType),'hermite')
    H = hermite_matrix(order);          % Hermite transformation matrix
    gain = 2*xstd/xrng;
    offset = 2*(xmean - xavg)/xrng;
    scale = [gain; offset];
    Ctemp = scale_poly(coeff,scale);
    Cnew = H\Ctemp;
else
    error('NewType not recognized');
end

  
Cnew = flipud(Cnew);



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














