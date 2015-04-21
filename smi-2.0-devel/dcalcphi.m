function [Phi,T] = dcalcphi(u,y,A,C,model)
 % dcalcphi  Construct the regressors matrix for the estimation of B, D 
%           and the initial state. optionally the Kalman gain regressors 
%           can be added when needed. This function is an internal 
%           function for the SMI toolbox. It is used in destbd, destx 
%           and dslslin. 
% 
% Syntax: 
%           [Phi,T]=dcalcphi(u,y,A,C,model) 
% 
% Input: 
%   u,y     The input and output data of the system to be identified. 
%   A,C     The system matrices A, C pair of the state space model 
%   model   five element flag vector [fB fD fx fK fY], indicating whether 
%           the regressors for B, D, x0 and K should be build. 
%           When fY is set, the stacked output is appended to the 
%           regressors matrix.  Default is [1 1 0 0 1]. 
% Output: 
%   Phi     The constructed regressors matrix. 
%   T       Transformation matrix of the internally used Schur form. 
%           After estimation use B=T*B, x0=T*x0, K=T*K and D=D'; 
% 
% See also: destbd, destx, dslslin 
 
%  --- This file is generated from the MWEB source destbd.web --- 
% 
% Bert Haverkamp, October 1999 
% copyright (c) 1999 Bert Haverkamp. 
 
 
 
 if nargin==0 
  help dcalcphi 
  return 
end 
if nargin<4 
  error('Not enough input variables.') 
end 
 
if size(y,2)>size(y,1) 
  y = y'; 
end 
if size(u,2)>size(u,1) 
  u = u'; 
end 
N = size(y,1); 
l = size(y,2); 
m = size(u,2); 
if ~(size(u,1)==N) 
  error('Input and output should have same lenght.') 
end 
if size(C,1)~=l 
  error('The number of outputs does not match the number of rows in C.') 
end 
 
error(abcdchk(A,[],C,[])); 
n = size(A,1); 
if (max(abs(eig(A)))>1) 
  disp('dcalcphi: A has instable pole. The estimate of B and D might be very bad.') 
end 
 
% test argument model 
if nargin<5 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0 0 1];% default 
end 
if length(model)~=5 
  error('Variable ''model'' should have five elements'); 
end 
 
fB = model(1); 
fD = model(2); 
fx = model(3); 
fK = model(4); 
fY = model(5); 
 
 
  
[T,A1] = schur(A); 
A1(abs(A1)<10 * eps) = 0;% arbitrary bound, this can break things! 
C1 = C * T; 
 
%Gamma 
if fx 
  temp = zeros(n * l,N); 
  Gamma = zeros(N * l,n); 
  Gamma(1:l,:) = C1; 
  An = A1; 
  for i = 1:floor(log(N)/log(2)), 
    Gamma(2^(i-1) * l+1:2^i * l,:) = Gamma(1:2^(i-1) * l,:) * An; 
    An = An * An; 
  end 
  Gamma(2^i * l+1:N * l,:) = Gamma(1:N * l-2^i * l,:) * An; 
  temp(:) = Gamma'; 
  for j = 1:l 
    Gamma(N * (j-1)+1:N * j,:) = temp(n * (j-1)+1:n * j,:)'; 
  end 
else 
  Gamma = []; 
end 
 
 
%Ybij 
if fB 
  Ybij = zeros(N * l,n * m); 
  e = eye(n); 
  for i = 1:n 
    if i==n 
      s = i; 
    elseif A1(i+1,i)==0 
      s = i; 
    else 
      s = i+1; 
    end 
    Aout = A1(1:s,1:s); 
    Bout = e(1:s,i); 
    Cout = C1(:,1:s); 
    for j = 1:m 
      x  =  ltitr(Aout,Bout,u(:,j),zeros(size(Aout,1),1)); 
      ybij  =  x  *  Cout.' ; 
      Ybij(:,(j-1) * n+i) = ybij(:); 
    end 
  end 
else 
  Ybij = []; 
end 
 
%Ykij 
if fK 
  Ykij = zeros(N * l,n * m); 
  e = eye(n); 
  for i = 1:n 
    if i==n 
      s = i; 
    elseif A1(i+1,i)==0 
      s = i; 
    else 
      s = i+1; 
    end 
    Aout = A1(1:s,1:s); 
    Bout = e(1:s,i); 
    Cout = C1(:,1:s); 
    for j = 1:l 
      x  =  ltitr(Aout,Bout,y(:,j),zeros(size(Aout,1),1)); 
      ykij  =  x  *  Cout.' ; 
      Ykij(:,(j-1) * n+i) = ykij(:); 
    end 
  end 
else 
  Ykij = []; 
end 
 
%Uij 
if fD 
  Uij = zeros(N * l,m * l); 
  for j = 1:l, 
    Uij((j-1) * N+1:j * N,(j-1) * m+1:j * m) = u; 
  end 
else 
  Uij = []; 
end 
 
%Y 
if fY 
  Y = zeros(l * N,1); 
  Y(:) = y; 
else 
  Y = []; 
end 
 
Phi = [Gamma,Ybij,Ykij,Uij,Y]; 
end




