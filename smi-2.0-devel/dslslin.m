 
function [A,B,C,D,x0,K,options]  =  dslslin(u,y,A,C,K,model,partype,options) 
 % dslslin  Performs a Least Squares optimization of a discrete 
%          time linear state space system system with model 
%          structure: 
% 
%              x(k+1) = Ax(k) + Bu(k) + Ke(k) 
%              y(k)   = Cx(k) + Du(k) +  e(k) 
% 
%          First, the state space matrices are parameterized. 
%          In order to minimize the number of parameters, the 
%          Separable Least Squares technique is used. This 
%          allows the cancellation of the parameters of B and D, 
%          such that only the parameters that describe A and 
%          C need to be optimized. The parameterized model is optimized 
%          with the leastsq function from the MATLAB 
%          optimization toolbox. If needed also the initial state 
%          and a Kalman gain can be optimized. 
% 
% Syntax: 
%          [A,B,C,D]=dslslin(u,y,A,C) 
%          [A,B,C,D,x0,K,options] = dslslin(u,y,A,C,K,model,partype,options) 
% Input: 
%  u,y     The input and output data of the system to be optimized. 
%  A,C     Initial estimates of the system matrices A, and C. 
%  model   Vector with flags that specify the kind of model to use. 
%          model(1) specifies if matrix B should be estimated. 
%          model(2) specifies if matrix D should be estimated. 
%          model(3) specifies if the initial state x0 should be estimated. 
%          model(4) specifies if the kalman filter gain K is estimated. 
%          Default is [1 1 0 0]. It is recommended only to 
%          estimate the initial state if the data length is 
%          short compared to the largest time constant of the system. 
% partype  This parameter specifies the type of parameterization 
%          that is used to parameterize the state space model. 
%          Two types of parameterization are supported: 
%          'on'= Output Normal and 'tr'=TRidiagonal. 
%  options Input parameters that are passed on directy to the 
%          optimization function from the Optimization Toolbox. 
%          See foptions for more information 
% 
% Output: 
%  A,B,C,D System matrices of the optimized linear model. 
%          If B or D is not estimated, it will be empty. 
%  x0      Estimate of the initial state. If the x0 matrix is 
%          not estimated, it will be returned empty. 
%  K       Estimate of the initial state. Same here. 
%  options Output parameters from the Optimization Toolbox. 
%          See foptions. 
% 
% See also:  leastsq, foptions 
 
%  --- This file is dslslin.m generated from the MWEB source dslslin.web --- 
% 
% B. Haverkamp and V. Verdult, August 1998 
% copyright (c) 1998 V. Verdult 
 
 
  
if nargin==0 
  help dslslin 
  return 
elseif nargin<4 
  error('Not enough input arguments') 
end 
 
%defaults: 
if nargin<8 
  options = []; 
end 
 
if nargin<7 
  partype = []; 
end 
if isempty(partype) 
  partype = 'on'; 
end 
 
if nargin<6 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0 0]; 
end 
 
if nargin<5 
  K = []; 
end 
 
 
if size(y,2)>size(y,1) 
  y = y'; 
end 
if size(u,2)>size(u,1) 
  u = u'; 
end 
 
m = size(u,2); 
l = size(y,2); 
N = size(u,1); 
params.N = N; 
params.m = m; 
params.l = l; 
 
if l==0, 
  error('We need an output') 
end 
if m==0, 
  error('We need an input') 
end 
 
error(abcdchk(A,[],C,[])) 
params.n = size(A,1); 
n = params.n; 
if size(C,1)~=l 
  Error('the number of outputs does not correspond with the size of C') 
end 
 
if length(model)~=4 
  error('The variable ''model'' must contain 4 values') 
else 
  fB = model(1); 
  fD = model(2); 
  fx = model(3); 
  fK = model(4); 
  params.fB = fB; 
  params.fD = fD; 
  params.fx = fx; 
  params.fK = fK; 
end 
 
if  partype=='on' 
  params.partype = partype; 
elseif  partype=='tr' 
  params.partype = partype; 
else 
  error('You specified an unknown type of parameterization in model') 
end 
 
if fK & isempty(K) 
  error('You specified fK=1, but did not provide an initial estimate'); 
end 
if (max(abs(eig(A)))>1) 
  error('Initial A must be stable.') 
end 
if rank(obsv(A,C))<n 
  error('Initial (A,C) must be observable.') 
end 
 
 
 
 
 
if fK 
  A = A-K * C; 
  % from now on, A is not really A. 
end 
thn = dss2th(A,C,params.partype); 
[thn,options] = leastsq('dfunlin',thn,options,[],u,y,params); 
[A,C] = dth2ss(thn,params); 
if isempty(A) 
  warning('Optimization did not result in stable model') 
  B = [];D = [];x0 = [];K = []; 
else 
  fB = params.fB; 
  fD = params.fD; 
  fx = params.fx; 
  fK = params.fK; 
  [B,D,x0,K] = destbd(u,y,A,C,[fB,fD,fx,fK]); 
  if fK 
    A = A+K * C; 
    % now A is the real A again 
    B = B+K * D; 
    else 
    K = x0; 
    x0 = options; 
  end 
end 
 
 

