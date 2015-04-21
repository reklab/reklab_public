 
function [A,B,C,D,x0,thl,options]  =  dslswie(u,z,A,B,C,D,x0,nn,model,  ... 
         partype,options) 
  
% dslswie  Performs a Least Squares optimization of a discrete 
%          time Wiener state space system system with model 
%          structure: 
%                 x(k+1) = Ax(k) + Bu(k) 
%                 y(k)   = Cx(k) + Du(k) 
%                 z(k)   = f(y(k)) 
%          First, the state space matrices are parameterized. 
%          In order to minimize the number of parameters, the 
%          Separable Least Squares technique is used. This 
%          allows the cancellation of the parameters of the nonlinearity, 
%          such that only the parameters that describe the linear part 
%          need to be optimized. The parameterized model is optimized 
%          with the leastsq function from the MATLAB optimization toolbox. 
% 
% Syntax: 
%          [A,B,C,D,x0,thl,options] = dslswie(u,z,A,B,C,D,x0,nn,model,options) 
% 
% Input: 
%  u,z     The input and output data of the system to be optimized. 
%  A,B,C,D Initial estimates of the linear part of the wiener model. 
%  x0      Initial condition, This is optional. 
%  nn      The order of the Chebychev polynomials in the static 
%          nonlinearity. 
%  model   Vector with flags that specify the kind of model to use. 
%          model(1) specifies if B should be estimated. 
%          model(2) specifies if D should be estimated. 
%          model(3) specifies if the initial state x0 should be 
%          estimated.  Default is [1 1 0]. It is recommended only 
%          to estimate the initial state if the data length is short 
%          compared to the largest time constant of the system. 
% partype  This parameter specifies the type of parameterization 
%          that is used to parameterize the state space model. 
%          Two types of parameterization are supported: 
%          'on'= Output Normal and 'tr'=TRidiagonal. 
% options  Input parameters that are passed on directly to the 
%          optimization function from the Optimization Toolbox. 
%          See foptions for more information. 
% 
% Output: 
% A,B,C,D  System matrices of the optimized linear model. 
%          If the D matrix is not estimated, it will be zero. 
% x0       Estimated initial state. If x0 is not estimated it \ 
%          will be empty. 
% options  Output parameters from the Optimization Toolbox. 
%          See foptions. 
% 
% See also:  foptions, dslslin 
 
%  --- This file is dslswie.m generated from the MWEB source dslswie.web --- 
% 
% B. Haverkamp and V. Verdult, August 1998 
% copyright (c) 1998 V. Verdult 
 
 
 
 if nargin==0 
  help dslswie 
  return 
elseif nargin<8 
  error('Not enough input arguments') 
end 
 
 
if nargin<11 
  options = []; 
end 
 
 
if nargin<10 
  partype = []; 
end 
if isempty(partype) 
  partype = 'on'; 
end 
if  partype=='on' 
  params.partype = partype; 
elseif  partype=='tr' 
  params.partype = partype; 
else 
  error('You specified an unknown type of parameterization in model') 
end 
 
 
% model = [fB, fD, fx] 
if nargin<9 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0]; 
end 
if length(model)~=3 
  error('model must contain 3 values') 
else 
  params.fB = model(1); 
  params.fD = model(2); 
  params.fx = model(3); 
end 
params.fK = 0;% no kalman filter in dslswie 
fB = params.fB; 
fD = params.fD; 
fx = params.fx; 
params.nn = nn; 
if fx+fB+fD==0 
  error('The system has no input or initial state. The output is zero') 
end 
 
if ~(max(abs(eig(A)))<1) 
 error('Initial A must be stable.') 
end 
if rank(obsv(A,C))<length(A) 
 error('Initial (A,C) must be observable.') 
end 
params.n = size(A,1); 
params.m = size(u,2); 
params.l = size(C,1); 
m = params.m; 
n = params.n; 
l = params.l; 
if fB & ~min(size(B)==[n,m]) 
  error('The size of B is incorrect'); 
end 
if fD & ~min(size(D)==[l,m]) 
  error('The size of D is incorrect'); 
end 
 
if size(z,2)>size(z,1) 
  z = z'; 
end 
if size(u,2)>size(u,1) 
  u = u'; 
end 
 
if size(z,2)==0, 
  error('We need an output') 
end 
if size(u,2)==0, 
  error('We need an input') 
end 
 
if ~(size(z,2)==l) 
  error('z must have as many columns as system outputs.') 
end 
N = size(z,1); 
if ~(size(u,1)==N) 
  error('Input and output must have the same length') 
end 
 
 
 thn = dss2th(A,B,C,D,x0,[],params.partype); 
[thn,options] = leastsq('dfunwie',thn,options,[],u,z,params); 
[A,B,C,D,x0] = dth2ss(thn,params); 
if ~fB 
  B = zeros(n,m); 
end 
if ~fD 
  D = zeros(l,m); 
end 
ye = dlsim(A,B,C,D,u); 
[thl,ze] = chebest(ye,z,nn); 
 
 
 
 

