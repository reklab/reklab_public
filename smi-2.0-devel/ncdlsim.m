function [y,x]=ncdlsim(Ac,Aa,Bc,Ba,Cc,Ca,D,u,xc0,xa0)
% ncdlsim   This function simulates a non-causal, unstable or unproper system.
%           The model is divided in a causal part, and an anti-causal 
%           part. The causal part is stable and can be simulated 
%           forward in time. The anti-causal part can only be simulated
%           stably backward in time.
%           The system can be given in the following forms
%  
%           Non causal state space form:
%               xc(k+1) = Ac xc(k) + Bc u(k)   (causal part)
%               xa(k-1) = Aa xa(k) + Ba u(k)   (anti-causal part)
%                y(k)   = Cc xa(k) + Ca xa(k) + D u(k)
%
%           Unstable state space form:
%                x(k+1) = A x(k) + B u(k) + w(k)
%                y(k)   = C x(k) + D u(k) + v(k)
%
%           Non-proper transfer function:
%           an*y(k+n)+ ... +a0*y(k) = bm*u(k+m)+ ... +b1*u(k+1)+ b0*u(k)
%
% Syntax:
%          [y,x]=ncdlsim(Ac,Aa,Bc,Ba,Cc,Ca,D,u,xc0,xa0)
%          [y,x]=ncdlsim(A,B,C,D,u,x0)
%          [y,x]=ncdlsim(num,den,u)
%
% Inputs:
% Ac,Aa,Bc,Ba State space matrices representing a non-causal state space 
% Cc,Ca,D     model.
% A,B,C,D     State space matrices representing a ordinary state space 
%             model, possibly unstable.
% num,den     Numerator and denominator of possibly non-proper
%             transfer function.
%             of possibly non-proper transfer function.
% u           Input to the  system
% xc0, xa0    Initial states of causal and anti-causal part of the
%             non-causal state space model.
% x0          Initial state of the ordinary state space model.
%
% Outputs:
%  y          Simulated output.
%  x          Simulated state of the state space model.
%             The first nc columns are the states of the causal part, 
%             the rest are the states of the anti-causal part.
%
% See also:  ncdestac, ncdestbd, dlsim

% Bert Haverkamp, 16-01-1995
% copyright(c)1995 Bert Haverkamp
% 16/07/96 BRJ Added transfer function support
% 09/10/96 BRJ Added unpropper transfer function support.
% 07/07/99 BRJ Rewrote everything, renamed to ncdlsim, added to smitoolbox

% for proper transfer function and state space we 
% first convert to the representation:
%                x(k+1) = A x(k) + B u(k) + w(k)
%                y(k)   = C x(k) + D u(k) + v(k)
%
% for the non proper transfer function we use 
% intermediate represtentation
%                x(k-1) = Ex(k) + B u(k)
%                y(k)   = Cx(k) + D(u(k)
if nargin==0
  help ncdlsim
  return
end


if (nargin==3)
  % transfer function
  num=Ac;den=Aa;u=Bc;
  rd=length(num)-length(den);% relative degree
  if rd>0 %non proper
    num=fliplr(num);%z->z^{-1}
    den=[zeros(1,rd),den];
    den=fliplr(den);
    [E,B,C,D]=tf2ss(num,den)
    proper_flag=0;
  else
    [A,B,C,D]=tf2ss(num,den);
    E=eye(size(A));
    x0=zeros(length(A),1);
    proper_flag=1;
  end

elseif (nargin==5)
  % state space model without initial state
  A=Ac;B=Aa;C=Bc;D=Ba;u=Cc;
  E=eye(size(A));
  x0=zeros(length(A),1);
  proper_flag=1;
elseif (nargin==6)
  % state space model with initial state
  A=Ac;B=Aa;C=Bc;D=Ba;u=Cc;
  E=eye(size(A));
  proper_flag=1;
end


if (nargin<8)
  if proper_flag
    % transform from intermediate to non-causal format
    [A,E,Q,Z,nc,na]=kroneckf(A,E);
    n=na+nc;
    Ac=A(1:nc,1:nc);
    Aa=E(nc+1:n,nc+1:n);
    B=Q*B;
    C=C*Z;
    Bc=B(1:nc,:);
    Ba=-B(nc+1:n,:);
    Cc=C(:,1:nc);
    Ca=C(:,nc+1:n)*Aa;
    if n>nc
      D=D-C(:,nc+1:n)*B(nc+1:n,:);
    end
    x0=inv(Z)*x0;
    xc0=x0(1:nc);
    xa0=zeros(n-nc,1); 
  else
    [A,E,Q,Z,nc,na]=kroneckf(eye(size(E)),E);
    n=na+nc;
    Ac=A(1:nc,1:nc);
    Aa=E(nc+1:n,nc+1:n);
    B=Q*B;
    C=C*Z;
    Bc=-B(1:nc,:);
    Ba=B(nc+1:n,:);
    Cc=C(:,1:nc)*Ac;
    Ca=C(:,nc+1:n);
    if nc>0
      D=D-C(:,1:nc)*B(1:nc,:);
    end
    xc0=zeros(nc,1);
    xa0=zeros(na,1);
  end
end

if nargin==8
  nc=length(Ac);
  na=length(Aa);
  xc0=zeros(nc,1);
  xa0=zeros(na,1);
end  


if size(u,2)>size(u,1)
  u=u';
end
m=size(u,2);

if isempty(Bc)
  Bc=zeros(0,m);
end
if isempty(Ba)
  Ba=zeros(0,m);
end
xc=ltitr(Ac,Bc,u,xc0);
xa=flipud(ltitr(Aa,Ba,flipud(u),xa0));
x=[xc,xa];
y=x*[Cc Ca]'+u*D';


