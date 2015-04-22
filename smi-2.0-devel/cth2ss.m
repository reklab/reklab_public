 function [A,B,C,D,x0,K]  =  cth2ss(theta,params,T) 
  
% cth2ss    This function converts a parameter vector 
%           that describes a continuous time state space model 
%           in output normal form to the state space matrices 
%           of that model. 
%           Model structure 
%              . 
%              x(t) = Ax(t) + Bu(t)+ K e(t) 
%              y(t) = Cx(t) + Du(t) + e(t) 
%              x(0) = x0 
% 
% Syntax: 
%           [A,B,C,D,x0,K] = cth2ss(theta,params) 
%           [A,C] = cth2ss(theta,params) 
% Input: 
%  theta    Parameter vector describing the system. 
%  params   A structure that contains the dimension parameters of 
%           the system, such as the order, the number of inputs, 
%           whether D, x0 or K are present, etc. 
% 
% Output: 
%  A,B,C,D  System matrices describing the state-space 
%           system in output normal form. If theta does 
%           not contain  parameters for D, this matrix 
%           will be returned as an empty matrix. 
%  x0       Initial condition. If theta does not contain 
%           parameters for x0 this vector will be returned 
%           as an empty matrix. 
%  K        Kalman gain. Same here. 
% 
% See also css2th, cslslin 
 
%  --- This file is generated from the MWEB source css2th.web --- 
% 
% Johan Bruls 1996, Bert Haverkamp 2000 
% copyright (c) 1996 Johan Bruls 
 
 
 
 if nargin==0 
  help cth2ss 
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
 
 
 
[A,C] = cthn2ac(thn,params); 
if nargout==2 
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
 function [A,C] = cthn2ac(thn,params); 
partype = params.partype; 
l = params.l; 
n = params.n; 
if partype=='on' 
  C  =  zeros(l,n);offset  =  0; 
  for j  =  1:l 
    C(j:l,j)  =  thn(offset+1:offset+l-j+1); 
    offset  =  offset+l-j+1; 
  end 
  Ass  =  zeros(n); 
  for j  =  1:min(l,n-1) 
    Ass  =  Ass+diag(thn(offset+1:offset+n-j),j)-diag(thn(offset+1:offset+n-j),-j); 
    offset  =  offset+n-j; 
  end 
  A  =  -.5  *  C'  *  C+Ass; 
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
  D(:) = thl(m * n+1:m * n+l * m); 
else 
  D = []; 
end 
if fx 
  x0 = thl(n * m+fD * m * l+1:n * m+fD * m * l+n); 
else 
  x0 = []; 
end 
if fK 
  K = zeros(n,l); 
  K(:) = thl(fx * n+m * n+fD * m * l+1:fx * n+m * n+fD * m * l+n * l); 
else 
  K = []; 
end 
 
  
 
 

