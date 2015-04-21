function [ze,Phi]=chebsim(y,thl)
% chebsim  Simulates a static nonlinear function on the basis of
%          Chebychev polynomials  with input y. The coeficients of
%          the Chebychev polynomials are given with the vector thl,
%          and commonly estimated with either dslswie or chebest.
%
% Syntax:
%          [ze,Phi]=chebsim(y,thl)
%
% Input:
% y        The input to the nonlinearity
% thl      Parameter vector, with coeficients of the nonlinear
%          function.
%
% Output:
% ze       Estimated output
% Phi      matrix with the Chebychev functions of y, such that 
%          ze=Phi*thl.
%
% see also: chebest, dslswie

% Bert Haverkamp, May 1996
% copyright (c) 1996 B.R.J. Haverkamp

if nargin==0
  help chebsim
  return
end
if nargin<2
  error('Not enough input variables')
end

if size(y,2)>size(y,1)
  y=y';
end

N=size(y,1);
l=size(y,2);

if size(thl,2)>size(thl,1)
  thl=thl';
end

if length(thl)<2*l
  error('Not enough parameters in thl')
end


% remove scaling and shifting coeficients first.
no_par=length(thl);
shift=thl(no_par-2*l+1:no_par-l);
scale=thl(no_par-l+1:no_par);
thl=thl(1:no_par-2*l);
no_lpar=length(thl);

if ~(round(no_lpar/l)==no_lpar/l)
  % does not detect every fault, just the obvious ones.
  error('wrong number of parameters in thl')
else
  nn=(no_lpar/l-1)/l;
end

if nn<1
  error('not enough parameters in thl')
end

Thl=zeros(1+nn*l,l);
Thl(:)=thl;
y=(y-ones(N,1)*shift')./(ones(N,1)*scale');
Phi=[ones(size(y)),y];
for j=3:nn+1
  Phi(:,(j-1)*l+1:j*l)=2*y.*Phi(:,(j-2)*l+1:(j-1)*l)-Phi(:,(j-3)*l+1:(j-2)*l);
end
Phi=Phi(:,l:(nn+1)*l);

ze=Phi*Thl; 





