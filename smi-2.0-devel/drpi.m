function [A,B,C,D]= drpi(u,y,i,n,beta,Ai,Bi,Ci,Di,Ri);
% drpi     Estimates in a recursive way the matrices 
%          of an LTI state space model in innovation form
%          using the output of the dordpi, destac and destbd
%          routines as initial guess. 
%          Model structure:
%              x(k+1) = Ax(k) + Bu(k) 
%              y(k)   = Cx(k) + Du(k) + v(k)
%          where v(k) is a zero-mean white noise sequence,
%          independent of the noise-free input u(k).
%                                       
% Syntax:
%          [A,B,C,D]= drpi(u,y,i,n,beta,Ai,Bi,Ci,Di,Ri);
%                                       
% Input:
%  u, y    The input and output data of the system to be identified.
%  i       The dimension parameter that determines the number of block
%          rows in the processed Hankel matrices such as ws used in dordpi
%          for the creation of Ri.
%  n       Order of system to be estimated. 
%  beta    Forgetting parameter for recursive algorithm (0 < beta < 1).
%  Ai,Bi,  Initial estimate of state space matrices.
%  Ci,Di   These matrices can for instance be obtained by using 
%          dordpi/destac/destbd on the first few samples.
%  Ri      R matrix from dordpi.
%                       
% Output:
% A,B,C,D  Estimated system matrices
%                                       
% See also: drpo, dordpi, destac, destbd

% The drpi routine is based on the PI scheme.

% Michel Verhaegen 1991
% Revised 1997 Marco Lovera


%
%       INITIALISATION
%
if nargin==0
  help drpi
  return
end

m=min(size(u));
mi = m*i;
mi2=2*mi;
l=min(size(y));
li=l*i;
li2=2*li;
mili=mi+li;
mi2li=mi2+li;
mi2li2=mi2+li2;

N=max(size(y));

beta2=beta^2;
betah=sqrt(beta);
lam2=beta;
J=[1 0;0 -1];

A=Ai;
C=Ci;
UU=[];

C2=100*eye(n,n);C2o=C2;

%
%       Initialisation of R
%

L0=Ri.L;
L=L0(:,1:2*mi);
NN=N-i+1-2*(m+l)*i;
L32=L((2*m)*i+1:2*m*i+l*i,m*i+1:(2*m)*i);
L32o=L32;
[U,S,Vt]=svd(L32);
Uo=U(:,1:n);

%
%       Initialisation of RLS for B and D
%

Th0=[reshape(Di,l*m,1);reshape(Bi,n*m,1)];

Thold=Th0;

P20=100*eye(l*m+n*m);

P2old=P20;

phi2t0=zeros(l,n*m);

for ii=1:i-2
  phi2t0=phi2t0+kron(u(ii,:),C*A^(i-ii-2));
  
end

phi2told=phi2t0;

%
%       MAIN LOOP
%

for s=2:NN,
  %
  %       UPDATE OF R
  %
  %
  %       Append column with new data
  %
  
  for ss=1:i,
    
    L((ss-1)*m+1:ss*m,mi2+1)=u(mili+s+i+ss-1,:)';   %Uf
    
    L((ss-1)*m+mi+1:ss*m+mi,mi2+1)=u(mili+s+ss-1,:)';       %Up
    
    L(mi2+(ss-1)*l+1:ss*l+mi2,mi2+1)=y(mili+s+i+ss-1,:)';   %Yf
    
  end
  
  %
  %       Givens rotations (part 1)
  %
  
  for ss=1:mi,
    modul=sqrt(L(ss,ss)^2 + L(ss,mi2+1)^2);
    co        = L(ss,ss)/modul;
    si        = L(ss,mi2+1)/modul;
    h         = co*L(ss:mi2li,ss) + si*L(ss:mi2li,mi2+1);
    L(ss:mi2li,mi2+1) = -si*L(ss:mi2li,ss) + co*L(ss:mi2li,mi2+1);
    L(ss:mi2li,ss)   = h;
  end
  
  phiy=L(mi2+1:mi2li,mi2+1);
  
  %
  %       Givens rotations (part 2)
  %
  
  for ss=mi+1:mi2,
    modul=sqrt(L(ss,ss)^2 + L(ss,mi2+1)^2);
    co        = L(ss,ss)/modul;
    si        = L(ss,mi2+1)/modul;
    h         = co*L(ss:mi2li,ss) + si*L(ss:mi2li,mi2+1);
    L(ss:mi2li,mi2+1) = -si*L(ss:mi2li,ss) + co*L(ss:mi2li,mi2+1);
    L(ss:mi2li,ss)   = h;
  end
  
  %
  %       Construction of phiybar 
  %
  
  phiybar=L(mi2+1:mi2li,mi2+1);
  
  %
  %       Update of C2 
  %
  
  X=[phiybar phiy];
  
  Y=Uo'*[phiybar phiy];
  
  temp=Y'*C2o;
  
  C2=1/beta*C2o-1/beta2*temp'*inv(1/beta*temp*Y+J)*temp;
  
  C2=(C2+C2')/2;
  
  %
  %       PA-Subspace tracking 
  %
  
  U=Uo+(X-Uo*Y)*J*Y'*C2;
  
  UU=[UU U];
  
  %
  %       Computation of A and C
  %
  
  U1   = U(1:(i-1)*l,1:n);
  
  U2   = U(l+1:i*l,1:n);
  
  A=[A U1\U2];
  
  C=[C U(1:l,1:n)]; 
  
  %
  %       Computation of B and D
  %
  
  %       Updating regressors for RLS
  
  phi2t=phi2told*kron(eye(m),A(:,(s-1)*n+1:s*n))+kron(u(s+i-2,:),C(:,(s-1)*n+1:s*n));
  
  phit=[kron(u(s+i-1,:),eye(l)) phi2t];
  
  %       RLS
  
  K2=P2old*phit'*inv(lam2*eye(l,l)+phit*P2old*phit');
  
  P2=1/lam2*(P2old-K2*phit*P2old);
  
  Th=Thold+K2*(y(s+i-1,:)'-phit*Thold);
  
  %       Extraction of B and D
  
  B(:,(s-1)*m+1:s*m)=reshape(Th(l*m+1:l*m+n*m),n,m);
  
  D(:,(s-1)*m+1:s*m)=reshape(Th(1:l*m),l,m);
  
  %
  %       Prepare the next iteration
  %
  
  Uo=U;
  
  for kk=1:2*mi
    L(kk:2*mi+li,kk)=L(kk:2*mi+li,kk)*betah;
  end
  C2o=C2;
  
  phi2told=phi2t;
  P2old=P2;
  Thold=Th;
  
end
