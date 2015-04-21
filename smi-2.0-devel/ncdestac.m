function [Ac,Aa,Cc,Ca]=ncdestac(R,n);
% ncdestac  Estimates the system matrices Ac, Aa, Cc and Ca of a
%           non-causal LTI state space model in innovation form 
%           using the output  of one of the dordxx preprocessor routines.
%           The model is divided in a causal part and an anti-causal 
%           part. The causal part is stable and can be simulated 
%           forward in time. The anti-causal part can only be
%           simulated stably backward in time.
%           The other matrices are estimated with ncdestbd. The model 
%           can be simulated with ncdlsim.
  
%           Model structure:
%               xc(k+1) = Ac xc(k) + Bc u(k)  (causal part)
%               xa(k-1) = Aa xa(k) + Ba u(k) (anti-causal part)
%                y(k)   = Cc xa(k) + Ca xa(k) + D u(k)
%
%           All preprocessor functions can be used for the causal 
%           case as well as the non-causal case. For information 
%           about the possible disturbance signals see the help pages 
%           of the preprocessor dordxx functions.
%
% Syntax:
% 	    [Ac,Aa,Cc,Ca]=ncdestac(R,n);
% 					
% Input:
%   R       Data structure obtained from dordxx, containing the 
%           triangular factor and additional information such as
%            (such as i/o dimension etc.)
%   n       Order of system to be estimated. 
% 			
% Output:
%  Ac,Aa,Cc,Ca  Estimated system matrices.
% 					
% See also: ncdlsim, ncdestbd 
if nargin==0
  help ncdestac
  return
end

m=R.m;l=R.l;i=R.i;L=R.L;Un=R.Un;
if (m<0)
  error('Illegal value for number of inputs in R matrix')
end
if (l<1)
  error('Illegal value for number of inputs in R matrix ')
end
if (i<0)
  error('Illegal value  for ''i'' in R matrix')
end

if ~(any(size(L,1)==[2*i*(m+l),2*i*(m+l)]))
%  error('L-matrix has unexpected size.')
end
if (n<1)
  error('System order of zero or lower does not make sense!')
end
if (n>(i-1)*l)
  error('n chosen too large, it should be smaller than or equal to ''(i-1)*l''.')
end

[U,S,V]=svd([Un(1:(i-1)*l,1:n) -Un(l+1:i*l,1:n)]);

A      = V(1:n,n+1:2*n);
E      = V(n+1:2*n,n+1:2*n);
[A,E,Q,Z,nc,na]=kroneckf(A,E);
U      = Un(:,1:n)*inv(Q);
Unc    = U(:,1:nc);
Una    = U(:,nc+1:n);

il = i*l;
if nc == 0    % anti-causal
  Ac=[];
  Cc=[];
  Ca = Una((i-1)*l+1:i*l,:);
  Aa  = Una(l+1:i*l,:)\Una(1:(i-1)*l,:);
elseif na==0, % causal
  Ac = Unc(1:il-l,1:n)\Unc(l+1:il,1:n);
  Cc = Unc(1:l,1:n);
  Aa=[];
  Ca=[];
else
  Ac = Unc(1:il-l,1:nc)\Unc(l+1:il,1:nc);
  Cc = Unc(1:l,1:nc);
  Aa = Una(l+1:il,1:na)\Una(1:il-l,1:na);
  Ca = Una(il-l+1:il,1:na);
end

if (max(abs(eig(Ac)))>1)|(max(abs(eig(Aa)))>1)
  % in some cases Ac and Aa can contain unstable poles
  % therefore we check and convert to stable Ac and Aa
  converted_flag=1;
  A=[Ac zeros(nc,na);zeros(na,nc),eye(na)];
  E=[eye(nc) zeros(nc,na);zeros(na,nc),Aa];
  C=[Cc,Ca];
  [A,E,Q,Z,nc,na]=kroneckf(A,E);
  Ac=A(1:nc,1:nc);
  Aa=E(nc+1:n,nc+1:n);
  C=C*Z;
  Cc=C(:,1:nc);
  Ca=C(:,nc+1:n)*Aa;
end





