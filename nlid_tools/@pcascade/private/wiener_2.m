function [irf, m,out] = wiener_2 (u,y,irflen,order,mode)
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

% Copyright 1992-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


data_len = length (u);


if nargin < 5
    mode = 'f'
  end

u = u(:);
y = y(:);
  
%  ignore the first irflen ponts, to avoid initial transients


y0 = y(irflen+1:data_len);
u0 = u(irflen+1:data_len);


%  use the principal generalized eigenvector of the second-order
%  input-output cross-correlation function and the Toeplitz structured input
%  autocorrelation matrix as an estimate of the impulse response of the
%  linear part of the Wiener system. 


phi_uu = phixy(u0,irflen);
PHIuu = toeplitz(phi_uu);
[U,S,V] = svd(PHIuu);
Sinv = diag(1./diag(S));
R = inv(chol(PHIuu));

%  Since the auto-correlation matrix is positive definite, use its Cholesky
%  factorization to reduce the problem to an ordinary eigenvalue problem


phi_uuy = phixxy (u0,y0,irflen);
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
[dummy,pi_order] = min(mdls);
irf = Hs(:,pi_order);



v = filter(irf,1,u);
v0 = v(irflen+1:data_len);
[m,vf] = tchebfit (v0,y0,order, mode);


%  If requested, calculate the output of the Wiener system due to the entire
%  input signal, initial transients and all.

if nargout > 2
  out = tchebval(m,v);
 end








