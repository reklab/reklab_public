function kern = gen_kern(basis,coeff,order);
%
% usage  kern = gen_kern(basis,coeff,order);
%
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


[hlen,N] = size(basis);

if order == 1
  kern = zeros(hlen,1);
  for i = 1:N
    kern = kern + coeff(i+1)*basis(:,i);
  end
elseif order == 2
  kern = zeros(hlen,hlen);
  coeffs = zeros(N,N);
  start = N+2;
  for i = 1:N
    stop = start + N - i;
    coeffs(i:N,i) = coeff(start:stop);
    coeffs(i,i:N) = coeff(start:stop)';
    start = stop + 1;
  end
  for v1 = 1:N
    kern = kern + coeffs(v1,v1)*basis(:,v1)*basis(:,v1)';
    for v2 = v1+1:N
      k1 = basis(:,v1); k2 = basis(:,v2);
      kern = kern +coeffs(v1,v2)*(k1*k2'+k2*k1')/2;
    end
  end
elseif order == 3
  kern = zeros(hlen^3,1);
  offset = (N+1)*(N+2)/2;
  for i = 1:N
    h1 = basis(:,i);
    for j = i:N
      h2 = basis(:,j);
      for k = j:N
	offset = offset+1;
	h3 = basis(:,k);
	kern = kern + coeff(offset)*k3_contrib(h1,h2,h3); 
      end
    end    
  end  
  kern = reshape(kern,hlen,hlen,hlen);
end


function k3 = k3_contrib(h1,h2,h3);

h1 = h1(:); h2 = h2(:); h3 = h3(:);
hlen = length(h1);
k3 =       reshape(reshape(h1*h2',hlen^2,1)*h3',hlen^3,1);
k3 = k3 +  reshape(reshape(h1*h3',hlen^2,1)*h2',hlen^3,1);
k3 = k3 +  reshape(reshape(h2*h3',hlen^2,1)*h1',hlen^3,1);
k3 = k3 +  reshape(reshape(h2*h1',hlen^2,1)*h3',hlen^3,1); 
k3 = k3 +  reshape(reshape(h3*h1',hlen^2,1)*h2',hlen^3,1);
k3 = k3 +  reshape(reshape(h3*h2',hlen^2,1)*h1',hlen^3,1);
k3 = k3/6;
   
