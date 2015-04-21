function y=tchebval(coeff,x,mode)
%
%	Function TCHEBVAL is used to evaluate the output
%	for a Tchebyschev polynomial with coefficients
%	COEFF, which can be obtained from the function
%	TCHEBFIT
%
%	USAGE : y=tchebval(coeff,x,mode)	 
%			coeff : chebyscev coefficients
%			x     : input
%			mode  : 'clip' or 'extend'
%
%	Note:  the vector coeff includes the domain minimum and 
%		maximum as the first two elements.  This is the 
%		only difference between this m-file and chebval
%
%		if only two arguments are passed, mode defaults 
%		to 'clip'.  The input is limited to the domain	
%		specified by the first two elements of coeff.
%
%
%  Calls: hard_limit,

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

if nargin == 2
    mode = 'clip';
    end
[nr,nc]=size(x);
x=x(:);
coeff = coeff(:);
len=length(coeff) - 2;
%
%  a and b are the minimum and maximum of the domain of definition of 
%  the polynomial.  The remainder of the vector "coeff" contains the 
%  coefficients of the tchebyshev polynomials.
%
a = coeff (1);
b = coeff (2);
coeff = coeff (3:len+2);
%
%  ifclipping has been enabled, clip the input to lie between a and b
%
if mode(1) == 'c'
    x = hard_limit (x,a,b);
    end
%
%  remap the input, such that a goes to -1 and b goes to +1
%
x=(2*x(:) - (b+a))/(b-a);
T=zeros(len,len);
A=zeros(length(x),len);
%
%  Create A, where the k'th column of A is the tchebyshev polynomial of 
%  order k-1, applied to the rescaled input vector x.  T contains the
%  tchebyshev coefficients.
%
T(1,1)=1;
T(2,1:2)=[1 0];
A(:,1)=polyval(T(1,1),x);
A(:,2)=polyval(T(2,1:2),x);
for i=3:len
	T(i,1:i)=[2*T(i-1,1:(i-1)) 0] - [0 0 T(i-2,1:(i-2))];
	A(:,i)=polyval(T(i,1:i),x);
end
if len == 1
  y = A(:,1)*coeff;
 else
  y=A*coeff;
 end
if nc > 1
	y=y';
end


