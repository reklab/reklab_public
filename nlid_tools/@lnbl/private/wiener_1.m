function [h,m,out]= wiener_1(z,hlen,order,mode)
%
% Syntax:  [h,m,vf,out]= wiener_1(u,y,hlen,order,mode)
%
%  

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

u = double(z(:,1));
y = double(z(:,2));

mdls = zeros(hlen,1);
data_length = length(u);

%  remove first hlen points to eliminate transients

% phiuu = phixy(u0,hlen);
% phiuy = phixy(u0,y0,hlen);

PhiUU=cor(z(:,1), 'Nlags',hlen);
PhiUY=cor(z,'Nlags',hlen); 
% phiuu=double(PhiUU);
% phiuy=double(PhiUY);
% 
% 
% 
% [U,S,V] = svd(toeplitz(phiuu));
% S = diag(S);
% Sinv = diag(1./S);
% 
% % automatic subspace dimension selection
% for i = 1:hlen
%   PHIuui = V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)';
%   h = PHIuui*phiuy;
%   testout = filter(h,1,u);
%   resid = y(hlen+1:data_length) - testout(hlen+1:data_length);
%   loss_func = std(resid)^2;
%   mdls(i)= loss_func * (1 + log(data_length-hlen)*i/(data_length-hlen));
% end
% 
% [val,pi_order] = min(mdls);
% PHIuui = V(:,1:pi_order)*Sinv(1:pi_order,1:pi_order)*U(:,1:pi_order)';
% h = PHIuui*phiuy;
% testout = filter(h,1,u);
% Z=cat(2,testout(hlen+1:data_length),y(hlen+1:data_length));

h=irf(z,'Nlags',hlen,'Mode',mode); 
yp = nlsim(h,z(:,1));
Z=cat(2,z(:,1),yp);

m=polynom(Z, 'OrderMax',order,'Type','tcheb');

% 
% [m,vf] = tchebfit(testout(hlen+1:data_length),y(hlen+1:data_length),order,mode);
% 
% if nargout > 2
%   out = tchebval(m,testout);
% end



