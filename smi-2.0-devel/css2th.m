 
function [theta,params,T]  =  css2th(varargin); 
 % css2th    This function converts a continuous time state space 
%           model to a parameter vector that describes the model. 
%           Model structure 
%              . 
%              x(t) = Ax(t) + Bu(t)+ K e(t) 
%              y(t) = Cx(t) + Du(t) + e(t) 
%              x(0) = x0 
% 
% Syntax: 
%           [theta,params,T] = css2th(A,B,C,D,x0,K,partype); 
%           [theta,params]  = css2th(A,B,C,D); 
%           [theta,params,T] = css2th(A,C); 
% 
% Input: 
% A,B,C,D   System matrices describing the state-space system. 
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
% theta     Parameter vector describing the system. 
% params    A structure that contains the dimension parameters of 
%           the system, such as the order, the number of inputs 
%           whether D, x0 or K are present, etc. 
% T         Transformation matrix between the input state space system 
%           and the state space system in the form described by theta 
%           (output normal or tridiagonal). 
% 
% See also : cth2ss, cslslin 
 
%  --- This file is generated from the MWEB source css2th.web --- 
% 
% Johan Bruls 1996 Bert Haverkamp 2000 
% copyright (c) 1996 Johan Bruls 
 
 
 
 if nargin==0 
  help css2th 
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
if partype=='on' 
  if max(real(eig(A)))>0 
    error('System must be stable.') 
  end 
  if max(real(eig(A)))>1-10 * eps 
    disp('A pole was found on imaginary axis.'); 
    error('I can''t reliably calculate the output normal parameterization.') 
  end 
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
 
 
 
 
 
[thn,T] = cac2thn(A,C,params); 
 
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
 function [thn,T] = cac2thn(A,C,params) 
partype = params.partype; 
n = params.n; 
l = params.l; 
if partype=='on' 
  % output normal 
  [A,C,T] = cac2on(A,C,params); 
  thn  =  []; 
  for j  =  1:min(n,l)% retrieve the columns of 
    thn  =  [thn;C(j:l,j)];% C, our gamma's 
  end 
  Ass  =  (A-A')/2;% calculate the skew-symmetric part of A 
   for j  =  1:min(l,n-1) 
    thn  =  [thn;diag(Ass,j)];% retrieve subdiagonals 
  end% containing alpha's 
elseif partype=='tr' 
  % tri-diagonal 
  [A,C,T] = cac2tr(A,C); 
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
 
 
 
 
 function [A,C,T] = cac2on(A,C,params) 
  n = params.n; 
  l = params.l; 
  % transform observability gramian into identity 
  Eo  =  lyap(A',C'  *  C); 
  T1  =  chol(Eo); 
  A  =  T1  *  A/T1; 
  C =  C/T1; 
  % transform C matrix into lower triangular form 
  [T2,R]  =  qr(C'); 
  A  =  T2'  *  A  *  T2; 
  C  =  R'; 
  % Transform A into bandsymetric matrix, (without disturbing the lower 
  % triangular property of C) 
  T3  =  eye(n); 
  if n>=l 
    for i  =  0:ceil((n-l)/l)-1 
      A21i  =  A(l  *  i+l+1:n,l  *  i+1:l  *  i+l); 
      [Q22i,R21i]  =  qr(A21i); 
      T3i  =  eye(n); 
      T3i(l  *  i+l+1:n,l  *  i+l+1:n)  =  Q22i; 
      T3  =   T3  *  T3i; 
      A  =  T3i'  *  A  *  T3i; 
    end 
    C  =  C  *  T3; 
  end 
  % transform gamma(1,1) and the alphas on the 
  % first superdiagonal to be positive. 
  if n>1, 
    T4  =  diag(cumprod(sign([C(1,1);diag((A-A')/2,1)]))); 
  else 
    T4  =  sign(C(1,1)); 
  end 
  A  =  T4  *  A  *  T4; 
  C  =  C  *  T4; 
  T = pinv(T1) * T2 * T3 * T4; 
 
 
   
 function [A,C,T] = cac2tr(A,C) 
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
 
 
 
 

