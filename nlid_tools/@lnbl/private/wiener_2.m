function h_est = wiener_2 (uy,hin)
%
%  usage: [irf, poly, out] = wiener2 (u,y,irflen,order,mode);
%
%  given input u, and output, y, as well as the memory length of the linear
%  part of the system, and the order of the static nonlinearity, this 
%  function returns a Wiener system based on the second order cross-correlation
%  between input and output.
%
%  mode is simply passed to the routine tchebfit.
%
%
% Calls: phixy, phixxy, tchebfit, tchebval

% Copyright 1992-2003, Robert E Kearney and  David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


global PINV_ORDER

irflen = get(hin,'nLags');
Ts = get(uy,'domainIncr');

u = double(uy(:,1));
data_len = length (u);
  
%  ignore the first irflen ponts, to avoid initial transients

uy0 = uy(irflen+1:data_len,:);
u0 = double(uy0(:,1));
y0 = double(uy0(:,2));

%  use the principal generalized eigenvector of the second-order
%  input-output cross-correlation function and the Toeplitz structured input
%  autocorrelation matrix as an estimate of the impulse response of the
%  linear part of the Wiener system. 


phi = cor;
set(phi,'kernOrder',1,'corType','covar','biasMode','biased','nLags',irflen);
phi_uu = double(nlident(phi,u0));

PHIuu = toeplitz(phi_uu);
[U,S,V] = svd(PHIuu);
Sinv = diag(1./diag(S));
R = inv(chol(PHIuu));

%  Since the auto-correlation matrix is positive definite, use its Cholesky
%  factorization to reduce the problem to an ordinary eigenvalue problem

set(phi,'kernOrder',2);
phi_uuy = double(nlident(phi,uy0));
RphiR = R'*phi_uuy*R;
[Vp,Dp] = eig(RphiR);
[lambda,pos] = max(abs(diag(Dp)));
h0 = R*Vp(:,pos);
phi_uy = sqrt(lambda/2)*PHIuu * h0;


%  project the irf onto the subspace spanned by the first i singular vectors
%  of PHIuu, evaluate the output of this new filter, and compute the MDL
%  cost funtion, where the number of free parameters is taken to be the
%  subspace dimension.

for i = 1:irflen
  Hs(:,i) = V(:,1:i)*Sinv(1:i,1:i)*U(:,1:i)'*phi_uy;
  testout = sign(lambda)*filter(Hs(:,i),1,u0).^2;
  testout = testout - mean(testout);
  loss_func = min(std(y0 - testout)^2,std(y0 + testout)^2);
  mdls(i)= loss_func * (1 + log(data_len)*i/data_len);
end
%plot(mdls);
%drawnow;
[dummy,PINV_ORDER] = min(mdls);
h_est = hin;
set(h_est,'dataSet', Hs(:,PINV_ORDER),'domainIncr',Ts);

