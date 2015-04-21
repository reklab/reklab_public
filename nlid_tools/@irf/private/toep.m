function x=toep(r,b)
% solve the Toeplitx matrix equation Tx=b.
% THIS FUNCTION WILL SOLVE THE MATRIX EQUATION Tx=b, WHERE
%  T IS A TOEPLITTZ MATRIX DESCRIBED BY THE VECTOR r.
%  THE VECTOR r  BE NORMALIZED SO THAT r(1) = 1.  
%  LEVINSON'S ALGORITHM IS USED.
%
% 	Usage : x=toep(r,b)

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

[rr,rc]=size(r);
[br,bc]=size(b);
n=length(r);
if rc > 1
   r = r';
end
if bc > 1
   b = b';
end
scale=r(1);
b=b/scale;
r=r/scale;
r=r(2:length(r));
y = zeros(size(b));
x = y;
y(1) = -r(1);
x(1) = b(1);
a = -r(1);
beta = 1;
for k=1:(n-1)
   beta = (1 - a^2)*beta;
   u=(b(k+1) - r(1:k)'*x(k:-1:1))/beta;
   x(1:k) = x(1:k) + u*y(k:-1:1);
   x(k+1) = u;
   if k  < n-1
	a = -(r(k+1) + r(1:k)'*y(k:-1:1))/beta;
        y(1:k) = y(1:k) + a*y(k:-1:1);
	y(k+1) = a;
   end
end

