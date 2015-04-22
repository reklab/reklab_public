function [B,D,x0,K,R,Phi] =  destbd(u,y,A,C,model,Rold,iv,niv); 
  
% destbd  Estimates the matrices B and D of the state space model 
% 
%              x(k+1) = A x(k) + B u(k) + w(k) 
%              y(k)   = C x(k) + D u(k) + v(k) 
%         using the knowledge of the pair A, C. This function is able 
%         to concatenate different input-output data batches, through 
%         the matrix R and Rold). B and D are calculated by solving a 
%         linear least squares problem. The function 
%         can handle errors-in-variables models by using a 
%         instrumental variables method. 
% Syntax: 
%         [B,D]=destbd(u,y,A,C); 
%         [B,D,x0,R,Phi]=destbd(u,y,A,C,[fB fD fx],Rold,iv,niv); 
% 
% Input: 
%   u, y  The input and output data of the system to be identified. 
%   A, C  The estimated system matrices A and C. 
%   model Four element flag vector [fB fD fx] indicating whether 
%         B, D and x0 should be estimated. The default value is [1 1 0]. 
%         The matrix B or D  can be assumed zero by setting 
%         fB or fD to zero. The calculation of x0 can be ommitted by 
%         setting fx to zero. However, x0 will not be assumed zero. 
%         It's influence will still be taken into account for the 
%         computation of B and D. However, not calculating x0 makes 
%         the calculation more robust. Variables that are not asked for 
%         will be returned as empty matrices. 
%   Rold  R matrix obtained from previous data batch. This variable 
%         can be used to process data in batches, or to combine data 
%         from different experiments. 
%   iv    Instrumental variable. For instance for the errors-in- 
%         variables case the iv can be the past output or for 
%         closed loop data, the reference input. 
%   niv   Three element vector describing the number of lagged iv's 
%         used. niv(1) is the number of lagged iv's, niv(2) the number 
%         of lagged past outputs, and niv(3) the number of lagged past 
%         inputs. If only one element is present, no past data is 
%         used. If not given or empty appropriate default values are used. 
%         When past data is used as instrumental variable, about half 
%         of the input- output data is used for the instrument. 
% 
% Output: 
%   B, D  The estimated system matrices B and D. 
%   x0    Estimated initial state of the system. 
%   R     Compressed data matrix, storing information on the calcuation of 
%         the matrices B and D in following destbd. Used when analyzing 
%         multiple input-output data sequences. 
%   Phi   Regressor matrix that was used for the Least squares estimate. 
% 
% See also: destac, destx 
 
%  --- This file is generated from the MWEB source destbd.web --- 
% 
% Bert Haverkamp, August 1999 
% copyright (c) 1999 Bert Haverkamp, Michel Verhaegen 
 
 
 
 if nargin==0 
  help destbd 
  return 
elseif nargin<4 
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
n = size(A,1); 
if ~(size(u,1)==N) 
  error('Input and output should have same lenght.') 
end 
 
% test argument model 
if nargin<5 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0 0];% default 
end 
if length(model)==3 
  model = [model,0]; 
end 
if length(model)~=4 
  error('Variable ''model'' should have four elements'); 
end 
fB = model(1); 
fD = model(2); 
fx = model(3); 
fK = model(4); 
fR = (nargout>fB+fx+fD+fK); 
fPhi = (nargout>fB+fx+fD+fK+fR); 
 
if nargin < 6, 
  Rold  =  []; 
end 
nparam = n+fB * n * m+fK * n * l+fD * m * l;% number of parameters 
nn = nparam+1;% size of Phi 

if ~(isempty(Rold) |all(size(Rold)==[nn,nn])) 
  error('R-matrix has unexpected size.') 
end 
 
if (N<nn+size(Rold,1)) 
  error('Not enough data points to find an estimate.') 
end 
 
 
if nargin<7 
  fIV = 0; 
else 
  fIV = 1; 
  if length(iv)~=size(iv,1) 
    iv = iv'; 
  end 
  [NN,p] = size(iv); 
end 
 
  % check if A unstable, if so change the algorithm 
if (max(abs(eig(A)))>1) 
    funstable = 1; 
    if fIV 
      % see Chou and Verhaegen 
      % subspace algorithms for the identification of 
      % multivariable dynamic errors-in-variables models. 
      % appendix B 
      [K,prec,mesg]  =  place(A',C',(0:n-1)/n); 
      K  =  K'; 
      A  =  A - K  *  C; 
      z  =  dlsim(A,K,C,zeros(l),y); 
      y  =  y - z; 
    end 
else 
  funstable = 0; 
end 
 
 
 
% Some explanation of variable names: 
% N_ are lengths of signals, N is the length of u and y 
% Niv is the length of the instrument iv. when past input/outputs are 
% used as instruments Nf is the length of the future part, Np is the 
% past part. 
% S_ is the number of columns in the intrumental variable matrix IV. 
% Siv is the number of columns for iv, Su for u, and Sy for y. 
 
if fIV 
  if nargin<8 
    % default values for instruments 
    Sy  =  ceil(4 * nparam/(2 * l+m)); 
    Su  =  ceil(Sy/2); 
    if isempty(iv), Siv  =  0; else, Siv  =  Su;  end 
  else 
    if length(niv)==1, niv = [niv(1),0, 0]; end 
    if length(niv)~=3, error('niv should contain 1 or 3 elements'); end 
    Siv = niv(1); 
    Sy = niv(2); 
    Su = niv(3); 
  end 
  if (Siv * p+Su * m+Sy * l)<nparam 
    disp('niv(1)*#iv +niv(2)*#outputs+niv(3)*#inputs should be bigger than') 
    disp(' the number of parameters that are estimated') 
    error('number of instrumental variables is chosen to small') 
 
  end 
  % compensate for hankel structure of IV 
  % such that following holds: Niv=N+Siv-1 
  % calculate size of past data (Np) and future data (Nf) 
  if max(Sy,Sy)>0,Nf1 = floor((N-max(Su,Sy))/2);else,Nf1 = N; end 
  if Siv>0,       Nf2 = N-Siv;                else,Nf2 = N; end 
  Nf = min(Nf1,Nf2); 
  Np = N-Nf; 
  up = u(1:Np,:);u = u(Np+1:N); 
  yp = y(1:Np,:);y = y(Np+1:N); 
end 
 
  
[Phi,T] = dcalcphi(u,y,A,C,[fB fD 1 fK 1]); 
% solve linear equation: Y=[Gamma Yij Uij] XBD 
% or something similar depending on the flags. 
 
% Construct IV 
if ~fIV 
  if ~isempty(Rold) 
    Rold = Rold(n+1:nn,:); % remove initial state part 
  end 
  R  =  triu(qr([Rold;Phi])); 
  R  =  R(1:nn,1:nn); 
else 
  IViv = zeros(Nf * l,p * Siv); 
  for i = 1:Siv 
    IViv(1:Nf,(i-1) * p+1:i * p) = iv(i:Nf+i-1,:); 
  end 
  %IVu 
  IVu = zeros(Nf * l,m * Su); 
  for i = 1:Su 
    IVu(1:Nf,(i-1) * m+1:i * m) = up(i:Nf+i-1,:); 
  end 
  % IVy 
  IVy = zeros(Nf * l,l * Sy); 
  for i = 1:Sy 
    IVy(1:Nf,(i-1) * l+1:i * l) = yp(i:Nf+i-1,:); 
  end 
  for j = 2:l 
    IViv((j-1) * Nf+1:j * Nf,:) = IViv(1:N,:); 
    IVu((j-1) * Nf+1:j * Nf,:) = IVu(1:Nf,:); 
    IVy((j-1) * Nf+1:j * Nf,:) = IVy(1:Nf,:); 
  end 
  IV = [IVy IVu IViv]; 
  % use small batch method from eiv_bd 
  Nk  =  min(round(500/nparam),Nf);  % #data points to be processed per batch 
  Nb  =  ceil(Nf/Nk);    % total number of batches 
  R = []; 
  for i = 1:Nb 
    Ns = (i-1) * Nk+1; 
    Ne = min(i * Nk,Nf); 
    Phi2 = IV(Ns:Ne,:)' * Phi(Ns:Ne,:); 
    R  =  triu(qr([Phi2;R])); 
    R  =  R(1:nn,1:nn); 
  end 
end 
 
sb = ~fx * n+1;% start of block taken from R 
eb = n+fB * m * n+fK * n * l+fD * m * l; % end of block from R 
yb = eb+1; 
XBD =  inv(R(sb:eb,sb:eb)) * R(sb:eb,yb); 
 
if fx 
  x0 = XBD(1:n); 
  x0 = T * x0; 
end 
if fB 
  B = zeros(n,m);B(:) = XBD(fx * n+1:fx * n+n * m); 
  B = T * B; 
end 
if fK 
  K = zeros(n,l);K(:) = XBD(fx * n+fB * n * m+1:fx * n+fB * n * m+n * l); 
  K = T * K; 
end 
if fD 
  D = zeros(m,l);D(:) = XBD(fx * n+fB * n * m+fK * n * l+1:fx * n+fB * n * m+fK * n * l+l * m); 
  D = D'; 
end 
 
if fK& fB &fD 
  B = B+K * D; 
end 
 
 
if funstable & fD &fB & fIV 
  B  =  B + K  *  D; 
end 
if ~fR 
  R = []; 
end 
if ~fK 
  K = []; 
end 
if ~fx 
  x0 = []; 
end 
if ~fD 
  D = []; 
end 
if ~fB 
  B = []; 
end 
 
 
 
if ~fK 
  K = R; 
  R = Phi; 
end 
 
% internal functions: 
 function [Phi,T] = dcalcphi(u,y,A,C,model); 
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
 
 
 
 
 
 
 
 

