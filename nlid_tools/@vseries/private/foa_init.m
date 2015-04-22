function [PmPr,G] = foa_init(u,y,hlen,Qmax);
%  performs initial computations for fast orthogonal algorithm, 
%  and subsequent fast basis expansions
%
%  syntax  [PmPr,G] = foa_init(u,y,hlen);
% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 
N = length(u);
M = prod([hlen+1:hlen+Qmax])/factorial(Qmax);
G = zeros(M,1);
%  M is the total number of regressors
%  mean and cross-correlation calculations
phi = cor;
set(phi,'kernOrder',1,'biasMode','biased','nSides',1,...
    'nLags',hlen,'corType','correl');
mu_y = mean(y);
G(1) = mu_y;
if Qmax > 0
  phiuu = double(nlident(phi,u));
  phiuy = double(nlident(phi,[u y]));
  G(2:hlen+1) = phiuy;
end
if Qmax > 1
  phi3u = double(nlident(phi,u,'kernOrder',2));
phi4u = double(nlident(phi,u,'kernOrder',3));
  phi4u = phi4u(:);
  phiuuy = double(nlident(phi,[u y],'kernOrder',2));
  
  m = hlen+2;
  for i = 1:hlen
    m1 = hlen - i;
    G(m:m+m1) = phiuuy(i,i:hlen);
    m = m + m1 + 1;
  end
end
switch Qmax
  case 0 
    PmPr = 1;
  case 1
    % we should really fix the pmpr mex file to do this,...
    % of course, if hlen is really big, we may run into memory problems.
    phi3u = zeros(hlen,hlen);
    phi4u = zeros(hlen^3,1);
    PmPr = fast_pmpr(phiuu,phi3u,phi4u,u);
    PmPr = PmPr(1:hlen+1,1:hlen+1);
  case 2
    PmPr = fast_pmpr(phiuu,phi3u,phi4u,u);
end
