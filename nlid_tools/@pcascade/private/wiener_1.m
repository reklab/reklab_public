function [h,m,out]= wiener_1(u,y,hlen,order,mode)
%
% Syntax:  [h,m,vf,out]= wiener_1(u,y,hlen,order,mode)
%
%  

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

u = u(:);
y = y(:);

mdls = zeros(hlen,1);
data_length = length(u);

%  remove first hlen points to eliminate transients
u0 = u(hlen+1:data_length);
y0 = y(hlen+1:data_length);
phiuu = phixy(u0,hlen);
phiuy = phixy(u0,y0,hlen);


[U,S,V] = svd(toeplitz(phiuu));
S = diag(S);
Sinv = diag(1./S);

% automatic subspace dimension selection
for i = 1:hlen
  PHIuui = V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)';
  h = PHIuui*phiuy;
  testout = filter(h,1,u);
  resid = y0 - testout(hlen+1:data_length);
  loss_func = std(resid)^2;
  mdls(i)= loss_func * (1 + log(data_length-hlen)*i/(data_length-hlen));
end

[val,pi_order] = min(mdls);
PHIuui = V(:,1:pi_order)*Sinv(1:pi_order,1:pi_order)*U(:,1:pi_order)';
h = PHIuui*phiuy;
testout = filter(h,1,u);

[m,vf] = tchebfit(testout(hlen+1:data_length),y0,order,mode);

if nargout > 2
  out = tchebval(m,testout);
end



