 
function [theta,params,T] = dss2th(varargin) 
  
% dss2th    This function converts a discrete time state space 
%           model to a parameter vector that describes the model 
%           Model structure: 
%                x(k+1) = Ax(k) + Bu(k) + Ke(k) 
%                y(k)   = Cx(k) + Du(k) + e(k) 
% Syntax: 
%           [theta,T,params] = dss2th(A,C,partype); 
%           [theta,T,params] = dss2th(A,B,C,D,partype); 
%           [theta,T,params] = dss2th(A,B,C,D,x0,K,partype); 
% 
% Input: 
% A,B,C,D   System matrices describing the state space system. 
%           The B and D matrices are optional and can be left out 
%           or given as an empty matrix to indicate it is not part 
%           of the model. 
% x0        Initial condition, This is optional. 
% K         Kalman gain. Also this matrix is optional. 
% partype   This parameter specifies the type of parameterization 
%           that is used to parameterize the state space model. 
%           Two types of parameterization are supported: 
%           'on'= Output Normal and 'tr'=TRidiagonal 
% 
% Output: 
% theta     Parameters vector describing the system. 
% params    A structure that contains the dimension of 
%           the system, such as the order, the number of inputs, 
%           whether D, x0 or K is present, etc. 
% T         Transformation matrix between the input state space 
%           system and the state space system in the form described 
%           by theta.  (the one that is constructed by {\dthtoss} 
% 
% See also: dth2ss 
 
%  --- This file is generated from the MWEB source css2th.web --- 
% 
% Johan Bruls 1996 Bert Haverkamp 2000 
% copyright (c) 1996 Johan Bruls 
 
 
 
 if nargin==0 
  help dss2th 
  return 
elseif nargin==1 
  error('There are not enough input arguments.') 
elseif nargin==2 
  fB = 0;fx = 0;fD = 0;fK = 0;m = 0; 
  A = varargin{1}; 
  C = varargin{2}; 
  partype = 'on'; 
elseif nargin==3, 
  fB = 0;fx = 0;fD = 0;fK = 0;m = 0; 
  A = varargin{1}; 
  C = varargin{2}; 
  partype = varargin{3}; 
elseif nargin==4 
  fK = 0;fx = 0; 
  partype = 'on'; 
  A = varargin{1}; 
  B = varargin{2}; 
  C = varargin{3}; 
  D = varargin{4}; 
  if isempty(B), fB = 0;else fB = 1;end 
  if isempty(D), fD = 0;else fD = 1;end 
  m = max(size(D,2),size(B,2)); 
elseif nargin==5 
  fK = 0;fx = 0; 
  A = varargin{1}; 
  B = varargin{2}; 
  C = varargin{3}; 
  D = varargin{4}; 
  partype = varargin{5}; 
  if isempty(B), fB = 0;else fB = 1;end 
  if isempty(D), fD = 0;else fD = 1;end 
  m = max(size(D,2),size(B,2)); 
elseif nargin==6 
  fK = 0; 
  A = varargin{1}; 
  B = varargin{2}; 
  C = varargin{3}; 
  D = varargin{4}; 
  x0 = varargin{5}; 
  partype = varargin{6}; 
  if isempty(B), fB = 0;else fB = 1;end 
  if isempty(D), fD = 0;else fD = 1;end 
  if isempty(x0), fx = 0;else fx = 1;end 
  m = max(size(D,2),size(B,2)); 
  partype = 'on'; 
elseif nargin==7 
  A = varargin{1}; 
  B = varargin{2}; 
  C = varargin{3}; 
  D = varargin{4}; 
  x0 = varargin{5}; 
  K = varargin{6}; 
  partype = varargin{7}; 
  if isempty(B), fB = 0;else fB = 1;end 
  if isempty(D), fD = 0;else fD = 1;end 
  if isempty(x0), fx = 0;else fx = 1;end 
  if isempty(K), fK = 0;else fK = 1;end 
  m = max(size(D,2),size(B,2)); 
else 
  error('There is a wrong number of input arguments') 
end 
[l,n] = size(C); 
 
if length(partype)>2 
  partype = partype(1:2); 
end 
if partype=='on' 
  nn = n * l; 
elseif partype=='tr' 
  nn = n * l+3 * n-2; 
else 
  error('You specified an unknown type of parameterization') 
end 
nl = fB * m * n+fD * m * l+fx * n+fK * n * l; 
params = struct('n',n,'m',m,'l',l,'fB',fB,'fD',fD,'fx',fx,'fK',fK,  ... 
       'partype',partype); 
if ~fB B = [];end; 
if ~fK K = [];end; 
if ~fx x0 = [];end; 
if ~fD D = [];end; 
 
if partype=='on' 
  if max(abs(eig(A)))>1 
    error('The system must be stable.') 
  end 
  if  rank(obsv(A,C))<length(A) 
    error('The system must be observable.') 
  end 
  if max(abs(eig(A)))>1-10 * eps 
    disp('A pole was found on the unit circle.'); 
    error('I can''t reliably calculate the output normal parameterization.') 
  end 
end 
 
 
[thn,T] = dac2thn(A,C,params); 
if nargin>3, 
  if params.fB, B = inv(T) * B;end 
  if params.fK, K = inv(T) * K;end 
  if params.fx, x0 = inv(T) * x0;end 
  thl = bdxk2thl(B,D,x0,K,params); 
else 
  thl = []; 
end 
theta = [thn;thl]; 
% internal functions: 
 function [thn,T] = dac2thn(A,C,params) 
partype = params.partype; 
n = params.n; 
l = params.l; 
Il = eye(l); 
if partype=='on' 
  % output normal 
  [A,C,T] = dac2on(A,C,params); 
  thn = zeros(n * l,1); 
  Z = [C;A]; 
  for i = 1:n 
    si = Z(n-i+2:n+l-i+1,n-i+1); 
    ri = Z(n-i+1,n-i+1); 
    ti = si' * si; 
    Ti = eye(n+l,n+l); 
    if ti>eps 
      Ti(n-i+1:n+l-i+1,n-i+1:n+l-i+1)  =  [-si,Il-(1-ri)/ti  *  si  *  si';ri,si']; 
    else 
      Ti(n-i+1:n+l-i+1,n-i+1:n+l-i+1)  =  [-si,Il;ri,transpose(si)]; 
    end 
    Z = Ti * Z; 
    thn(l * (i-1)+1:l * i) = si; 
  end 
elseif partype=='tr' 
  % tri-diagonal 
  [A,C,T] = dac2tr(A,C); 
  thn = [diag(A,1);diag(A);diag(A,-1);C(:)]; 
end 
 
 
 
 function thl = bdxk2thl(B,D,x0,K,params) 
fD = params.fD; 
fx = params.fx; 
fK = params.fK; 
thl = B(:); 
if fD, thl = [thl;D(:)]; end 
if fx, thl = [thl;x0];   end 
if fK, thl = [thl;K(:)]; end 
 
 
 
 
 function [A,C,T] = dac2on(A,C,params) 
n = params.n; 
l = params.l; 
% Make system output normal 
Wo = dgram(A',C'); 
Q1  = chol(Wo); 
A = Q1 * A/Q1; 
C = C/Q1; 
Z = [C;A]; 
% Make observability matrix positive lower triangular 
Gam = obsv(A,C); 
[Q2,R] = qr(Gam'); 
%Gam*Q2 
Q3 = diag(sign(diag(R(1:n,1:n)))); 
Q4 = Q2 * Q3; 
A = Q4 \ A * Q4; 
C = C * Q4; 
T = pinv(Q1) * Q4; 
 
 
 
 function [A,C,T] = dac2tr(A,C) 
[V,Lambda]  =  eig(A); 
Lambda  =  diag(Lambda); 
% First, separate the real and complex eigenvalues. 
n  =  length(Lambda); 
flag  =  zeros(n,1); 
for kl  =  1:n 
  flag(kl)  =  isreal(Lambda(kl)); 
end 
RealIndices  =  find(flag); 
ComplexIndices  =  find(1-flag); 
nc  =  length(ComplexIndices); 
nr  =  length(RealIndices); 
LambdaReal  =  Lambda(RealIndices); 
LambdaComplex  =  Lambda(ComplexIndices); 
Lambda  =  [LambdaComplex; LambdaReal]; 
VReal  =  V(:,RealIndices); 
VComplex  =  V(:,ComplexIndices); 
V  =  [VComplex VReal]; 
% Now, convert pairs of eigenvalues into 2by2 blocks -- except for the last 
% (real) eigenvalue, when n is odd. 
T  =  eye(n); 
for kl  =  1:2:nc 
  lambda  =  Lambda(kl); 
  a  =  real(lambda); 
  b  =  imag(lambda); 
  Mk  =  [a b;-b a]; 
  [Vk,Dk]  =  eig(Mk); 
  Tk  =  inv(Vk); 
  T(kl:kl+1,kl:kl+1)  =  Tk; 
end 
for kl = nc+1:2:nc+nr-1 
  lambda1  =  Lambda(kl); 
  lambda2  =  Lambda(kl+1); 
  a  =  (lambda1 + lambda2)/2; 
  b  =  (lambda1 - lambda2)/2; 
  Mk  =  [a b; b a]; 
  [Vk,Dk]  =  eig(Mk); 
  Tk  =  inv(Vk); 
  T(kl:kl+1,kl:kl+1)  =  Tk; 
end 
T  =  V * T; 
A  =  inv(T) * A * T; 
A  =  tril(triu(A,-1),1); 
C  =  C * T; 
 
 
 
 
 

