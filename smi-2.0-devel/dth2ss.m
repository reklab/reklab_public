 function [A,B,C,D,x0,K] = dth2ss(theta,params,T) 
  
% dth2ss   This function converts a parameter vector 
%          that describes a discrete time state space model 
%          in output normal form to the state space matrice 
%          of that model. 
%          Model structure: 
%               x(k+1) = Ax(k) + Bu(k) + Ke(k) 
%               y(k)   = Cx(k) + Du(k) + e(k) 
% Syntax: 
%          [A,B,C,D,x0,K] = dth2ss(theta,params) 
%          [A,C] = dth2ss(theta,params) 
% Input: 
% theta    Parameter vector describing the system. 
% params   A structure that contains the dimension parameters of 
%          the system, such as the order, the number of inputs, 
%          whether D, x0 or K is present in the model, etc. 
% T        Transformation matrix to be applied to the state space 
%          system that is constructed from theta. This transformation 
%          might come from the function dss2th and can be used to 
%          reconstruct the original state space matrices that were given 
%          to dss2th. 
% 
% Output: 
% A,B,C,D  System matrices describing the state space 
%          system in output normal form. If theta does 
%          not contain  parameters for  a matrix, this matrix 
%          will be returned as an empty matrix. 
%   x0     Initial state. If theta does not contain 
%          parameters for x0, this vector will be returned 
%          as an empty matrix. 
%   K      Kalman gain. Same here. 
% 
% See also: dss2th 
 
%  --- This file is generated from the MWEB source dss2th.web --- 
% 
% Johan Bruls 1996, Bert Haverkamp 2000 
% copyright (c) 1996 Johan Bruls 
 
 
 
 if nargin==0 
  help dth2ss 
  return 
elseif nargin==1 
  error('There are not enough input arguments.') 
end 
 
%defaults 
if nargin<3 
  T = []; 
end 
 
if isstruct(params) 
  n = params.n; 
  l = params.l; 
  m = params.m; 
  if nargout==2 
    params.fB = 0; 
    params.fD = 0; 
    params.fx = 0; 
    params.fK = 0; 
  end 
  fB = params.fB; 
  fD = params.fD; 
  fx = params.fx; 
  fK = params.fK; 
  partype = params.partype; 
end 
 
if partype=='on' 
  nn = n * l; 
elseif partype=='tr' 
  nn = n * l+3 * n-2; 
else 
  error('You specified an unknown type of parameterization in params.partype') 
end 
nt = length(theta); 
nl = fx * n+fB * n * m+fD * l * m+fK * n * l; 
if (nt~=nn+nl) 
    error(['The length of theta and the values in params do not'  ... 
    ' correspond.']); 
end 
thn = theta(1:nn); 
thl = theta(nn+1:nt); 
 
 
[A,C] = dthn2ac(thn,params); 
if nargout<=2 
  B = C; 
else 
  [B,D,x0,K] = thl2bdxk(thl,params); 
end 
 
if ~isempty(T) 
  A = T * A * inv(T); 
  C = C * inv(T); 
  if fB, B = T * B;   end 
  if fx, x0 = T * x0; end 
  if fK, K = T * K;   end 
end 
 
% internal functions: 
 function [A,C] = dthn2ac(thn,params) 
partype = params.partype; 
n = params.n; 
l = params.l; 
if partype=='on' 
  % output normal 
  Il = eye(l); 
  Z = [zeros(l,n);eye(n)]; 
  for i = 1:n 
    si = thn(l * (n-i)+1:l * (n-i+1)); 
    ti = si' * si; 
    if ti>1 
      % in case of invalid vector, return empty A matrix and return 
      warning('invalid theta vector') 
      A = []; 
      C = []; 
      return 
    end 
    ri = sqrt(1-ti); 
    Ti = eye(n+l,n+l); 
    if ti>eps 
      Ti(i:l+i,i:l+i)  =  [-si,Il-(1-ri)/ti * si * si';ri,si']; 
    else 
      %limit ti->0 
      Ti(i:l+i,i:l+i)  =  [-si,Il-0.5 * si * si';ri,si']; 
    end 
    Z = Ti' * Z; 
  end 
  C = Z(1:l,:); 
  A = Z(l+1:l+n,:); 
elseif partype=='tr'; 
  % tri-diagonal 
  A = zeros(n); 
  C = zeros(l,n); 
  for j = 1:n-1 
    A(j,j+1) = thn(j); 
    A(j,j) = thn(n-1+j); 
    A(j+1,j) = thn(2 * n-1+j); 
  end 
  A(n,n) = thn(2 * n-1); 
  C(:) = thn(3 * n-1:3 * n-2+n * l); 
end 
 
 
 function [B,D,x0,K] = thl2bdxk(thl,params) 
fB = params.fB; 
fx = params.fx; 
fD = params.fD; 
fK = params.fK; 
n = params.n; 
m = params.m; 
l = params.l; 
 
 
if fB 
  B = zeros(n,m); 
  B(:) = thl(1:m * n); 
else 
  B = []; 
end 
 
if fD 
  D = zeros(l,m); 
  D(:) = thl(fB * m * n+1:fB * m * n+l * m); 
else 
  D = []; 
end 
if fx 
  x0 = thl(fB * n * m+fD * m * l+1:fB * n * m+fD * m * l+n); 
else 
  x0 = []; 
end 
if fK 
  K = zeros(n,l); 
  K(:) = thl(fx * n+fB * m * n+fD * m * l+1:fx * n+fB * m * n+fD * m * l+n * l); 
else 
  K = []; 
end 
 
  
 
 

