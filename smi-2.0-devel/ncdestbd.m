function [Bc,Ba,D,xc0,xa0,R,Phi]= ncdestbd(u,y,Ac,Aa,Cc,Ca,model,Rold,iv,niv);
% ncdestbd  Estimates the matrices Bc, Ba and D of a non-causal 
%           state space model, using the knowledge of the matrices
%           Ac, Aa, Cc and Ca. This function is able  to concatenate 
%           different data batches, through the matrix R (R, Rold).
%           Model structure:
%               xc(k+1) = Ac xc(k) + Bc u(k)  (causal part)
%               xa(k-1) = Aa xa(k) + Ba u(k) (anti-causal part)
%                y(k)   = Cc xa(k) + Ca xa(k) + D u(k)
% 
%           Bc, Ba and D  are calculated by solving a linear least squares
%           problem. The influence of the initial state x0 in the data
%           batches is compensated for in the solution.  x0c
%           and x0a can also be estimated.
%
% Syntax:
%          [Bc,Ba,D,x0c,x0a,R]=ncdestbd(u,y,Ac,Aa,Cc,Ca,model,Rold,iv,niv);
%          [Bc,Ba,D]=ncdestbd(u,y,Ac,Aa,Cc,Ca);
% Input:
%  u,y      The input respectively output data of the system to be identified.
%  Ac,Ae    The estimated system matrices Ac, Ae, Cc and Ca of the
%  Cc,Ca    non-causal state space system matrices.  
%  model    Two element vector [fB fD fx] indicating whether Bc, Ba, D
%           and x0c, x0a should be estimated. Default value is [1 1 0].
%  Rold     R matrix obtained from previous data batch. This variable 
%           can be used to process data in batches, or to combine data 
%           from different experiments.
%  iv       instrumental variable, for instance for the errors-in-variables
%           case the iv can be the past output or for closed loop data,
%           the reference input. 
%  niv      Three element vector describing the number of lagged iv's 
%           used. niv(1) is the number of lagged iv's, niv(2) the number 
%           of lagged past outputs, and niv(3) the number of lagged past 
%           inputs. If only one element is present, no past data is 
%           used. If not given or empty appropriate default values are used. 
%           When past data is used as instrumental variable, about half 
%           of the input- output data is used for the instrument. 

%
% Output:
%  Bc,Ba,D  The estimated system matrices Bc, Ba and D.
%  x0c,x0c  The estimated initial state of the system.
%  R        Compressed data matrix, storing information on the
%           calculation of the matrices Bc, Ba and D in a following ncdestbd.
%	    Used when analyzing multiple input-output data sequences.
%
% See also:  ncdlsim, ncdestac


% Bert Haverkamp, Oct 2000
% copyright (c) 2000 Bert Haverkamp
  
if nargin==0
  help ncdestbd
  return
end

if nargin<6
  error('Not enough input variables.')
end

if size(y,2)>size(y,1)
  y=y';
end
if size(u,2)>size(u,1)
  u=u';
end
nc=size(Ac,1);
na=size(Aa,1);
n=nc+na;
N=size(y,1);
l=size(y,2);
m=size(u,2);
if nc>0
  if size(Cc,1)~=l
    error('The number of outputs and number of rows in Cc don''t match.')
  end
end
if na>0
  if size(Ca,1)~=l
    error('The number of outputs and number of rows in Ca don''t match.')
  end
end

if ~(size(u,1)==N)
  error('Input and output should have same lenght.')
end

if (max(abs(eig(Ac)))>1) |(max(abs(eig(Aa)))>1)
  funstable=1;
else
  funstable=0;
end

%test argument model 
if nargin<7
  model = []; 
end 

% test argument model 
if isempty(model) 
  model = [1 1 0];% default 
end 

if length(model)~=3 
  error('Variable ''model'' should have three elements'); 
end 
fB = model(1); 
fD = model(2); 
fx = model(3); 
fR = (nargout>2*fB+2*fx+fD); 
fPhi = (nargout>2*fB+2*fx+fD+fR); 

if nargin < 8, 
  Rold = []; 
end


nparam = n+fB * n * m+fD * m * l;% number of parameters 
nn = nparam+1;% size of Phi 


if ~(isempty(Rold) |all(size(Rold)==[nn,nn])) 
  error('R-matrix has unexpected size.') 
end 

if (N<nn+size(Rold,1)) 
  error,('Not enough data points to find an estimate.') 
end 
 
if nargin<9
  fIV = 0; 
else 
  fIV = 1; 
  if length(iv)~=size(iv,1) 
    iv = iv'; 
  end 
  [NN,p] = size(iv); 
end 



% Some explanation of variable names: 
% N is lengths of signals,
% Niv is the length of the instrument iv. 
% When past input/outputs are  used as instruments N is changed 
% to the length of the future part, Np is the past part. 
% Siv is the number of columns for iv, Su for u, and Sy for y. 
if fIV 
  if nargin<10
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
  if Siv>0,       Nf2 = N-Siv;                  else,Nf2 = N; end 
  Nf = min(Nf1,Nf2); 
  Np = N-Nf; 
  up = u(1:Np,:);u = u(Np+1:N); 
  yp = y(1:Np,:);y = y(Np+1:N); 
  N=Nf;
end 
 
if funstable
  % convert to stable Ac and Aa
  A=[Ac zeros(nc,na);zeros(na,nc),eye(na)];
  E=[eye(nc) zeros(nc,na);zeros(na,nc),Aa];
  C=[Cc,Ca];
  [A,E,Q,Z,nc,na]=kronekf(A,E);
  Ac=A(1:nc,1:nc);
  Aa=E(nc+1:n,nc+1:n);
  C=C*Z;
  Cc=C(:,1:nc);
  Ca=C(:,nc+1:n)*Aa;
end

% Y
Y=zeros(l*N,1);
Y(:)=y;

[Qc,Ac1]=schur(Ac);
Cc1=Cc*Qc;
[Qa,Aa1]=schur(Aa);
Ca1=Ca*Qa;


%Gammac
if nc>0
  temp=zeros(nc*l,N);
  Gammac=zeros(N*l,nc);
  Gammac(1:l,:)=Cc1;
  An=Ac1;
  for i=1:floor(log(N)/log(2)),
    Gammac(2^(i-1)*l+1:2^i*l,:)=Gammac(1:2^(i-1)*l,:)*An;
    An=An*An;
  end
  Gammac(2^i*l+1:N*l,:)=Gammac(1:N*l-2^i*l,:)*An;
  temp(:)=Gammac';
  for j=1:l
    Gammac(N*(j-1)+1:N*j,:)=temp(nc*(j-1)+1:nc*j,:)';
  end
else
  Gammac=[];
end


%Gammaa
if na>0
  temp=zeros(na*l,N);
  Gammaa=zeros(l*N,na);
  Gammaa(l*N-l+1:l*N,:)=Ca1;
  An=Aa1;
  for i=1:floor(log(N)/log(2)),
    Gammaa(N*l-2^i*l+1:N*l-2^(i-1)*l,:)= Gammaa(N*l-2^(i-1)*l+1:N*l,:)* An;
    An=An*An;
  end
  Gammaa(1:N*l-2^i*l,:)=Gammaa(2^i*l+1:N*l,:)*An;
  temp(:)=Gammaa';
  for j=1:l
    Gammaa(N*(j-1)+1:N*j,:)=temp(na*(j-1)+1:na*j,:)';
  end
else
  Gammaa=[];
end  

%Uij
if fD
  Uij=zeros(N*l,m*l);
  for j=1:l,
    Uij((j-1)*N+1:j*N,(j-1)*m+1:j*m)=u;
  end
else
  Uij=[];
end

%Ycij
if fB
  Ycij=zeros(N*l,nc*m);
  e=eye(nc);
  temp = zeros(N*l,1);
  for i=1:nc 
    if i==nc
      s=i;
    elseif Ac1(i+1,i)==0
      s=i;
    else
      s=i+1;
    end
    Aout=Ac1(1:s,1:s);
    Bout=e(1:s,i);
    Cout=Cc1(:,1:s);
    for j=1:m
      x = ltitr(Aout,Bout,u(:,j),zeros(size(Aout,1),1));
      ycij = x * Cout.' ;  
      Ycij(:,(j-1)*nc+i)=ycij(:);
    end
  end
else
  Ycij=[];
end

%Yaij
if fB
  Yaij=zeros(N*l,na*m);
  e=eye(na);
  temp = zeros(N*l,1);
  for i=1:na
    if i==na
      s=i;
    elseif Aa1(i+1,i)==0
      s=i;
    else
      s=i+1;
    end
    Aout=Aa1(1:s,1:s);
    Bout=e(1:s,i);
    Cout=Ca1(:,1:s);
    for j=1:m
      x=flipud(ltitr(Aout,Bout,flipud(u(:,j)),zeros(size(Aout,1),1)));
      yaij = x * Cout' ;  
      temp(:)=yaij;
      Yaij(:,(j-1)*na+i)=yaij(:);
    end
  end
else
  Yaij=[];
end

%IV 
if fIV 
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
end 

% resulting linear equation: 
% Y =[Gammac Gammaa Ycij Yaij Uij]*XBD

Phi=[Gammac,Gammaa,Ycij,Yaij,Uij,Y];

if fIV 
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
else
  if ~isempty(Rold)
    Rold=Rold(n+1:nn,:);  % remove initial state part
  end
  R = triu(qr([Rold;Phi,Y]));
  R = R(1:nn,1:nn);
end
sb = ~fx * n+1;% start of block taken from R 
eb = n+fB * m * n+fD * m * l; % end of block from R 
yb = eb+1; 
XBD =  inv(R(sb:eb,sb:eb)) * R(sb:eb,yb); 

if fx 
  xc0=XBD(1:nc);
  xc0=Qc*xc0;
  xa0=XBD(nc+1:n);
  xa0=Qa*xa0;
  if funstable
    x0=inv(Q)*[xc0;xa0];
    xc0=x0(1:nc);
    xa0=x0(nc+1:n);
  end
end

if fB 
  Bc=zeros(nc,m);Bc(:)=XBD(fx*n+1:fx*n+nc*m);
  Bc=Qc*Bc;
  Ba=zeros(na,m);Ba(:)=XBD(fx*n+nc*m+1:fx*n+n*m);
  Ba=Qa*Ba;	
  if funstable
    B=inv(Q)*[Bc;Ba];
    Bc=B(1:nc,:);
    Ba=B(nc+1:n,:);
  end
end 
if fD 
  D=zeros(m,l);D(:)=XBD(fx*n+fB*n*m+1:fx*n+fB*n*m+l*m);
  D=D';
end 


if ~fR 
  R = []; 
end 
if ~fx 
  xc0 = []; 
  xa0 = []; 
end 
if ~fD 
  D = []; 
end 
if ~fB 
  Bc=[]; 
  Ba=[];
end 
 