function [thl,ze,Phi]=chebest(y,z,nn)
% chebest   This function estimates a MIMO static nonlinear function
%           between the signals y and z. The function is
%           estimated on the basis of Chebychev polynomials.
%           Before estimating the coefficients of the polynomials
%           the input signal is shifted and scaled to fall within the
%           region [-1,1]. The shifting and scaling factors are
%           included in the parameter vector.
%
% Syntax:
% 	    [thl,ze,Phi]=chebest(y,z,nn)
% 	   
% Input:
%   y,z     Input and output of the nonlinearity.
%   nn      Order of the Chebychev polynomials in the nonlinear
%           function.
%                                       
% Output:
%   thl     Vector with the parameters of the static nonlinearity.
%           the last two elements of thl are the shift and scaling  
%           that was applied to y.
%   ze      Estimated output, on basis of the model that is obtained.
%   Phi     matrix with the Chebychev functions of y, such that 
%           ze=Phi*thl.
  
%
% see also: chebsim, dslswie


% Bert Haverkamp, May 1996
% copyright (c) 1996 B.R.J. Haverkamp

if nargin==0
  help chebest
  return
end
if nargin<3
  error('Not enough input variables')
end

if size(y,2)>size(y,1)
  y=y';
end

if size(z,2)>size(z,1)
  z=z';
end

N=size(y,1);
l=size(y,2);


if ~(size(z,1)==N)
  error('Both signals should have same lenght')
end

if ~(size(z,2)==l)
  error('Both signals should have same number of elements')
end

if nn<1, 
  error('nn should be positive')
end


shift=(max(y)+min(y))/2;
scale=(max(y)-min(y))/2;
y=(y-ones(N,1)*shift)./(ones(N,1)*scale);
Phi=[ones(size(y)),y];% zero and first order 
for j=3:nn+1
  Phi(:,(j-1)*l+1:j*l)=2*y.*Phi(:,(j-2)*l+1:(j-1)*l)-Phi(:,(j-3)*l+1:(j-2)*l);
end
Phi=Phi(:,l:(nn+1)*l);% improve condition number by removing first l-1 columns

thl=pinv(Phi)*z;
ze=Phi*thl; 

if nargout==0 % plot the nonlinearity
  figure
  for k=1:l;
    %  x is a linear line that spans the whole input range.
    x1=min(y(:,k)); x2=max(y(:,k)); step=(x2-x1)/N;
    x=zeros(N,l); x(:,k)=(x1:step:x2-step)';
    Phidum=[ones(size(x)),x];
    for j=3:nn+1
       Phidum(:,(j-1)*l+1:j*l)=2*x.*Phidum(:,(j-2)*l+1:(j-1)*l)-Phidum(:,(j-3)*l+1:(j-2)*l);
    end
    Phidum=Phidum(:,l:(nn+1)*l);

    zdum=Phidum*thl;
    for i=1:l
      subplot(l,l,(i-1)*l+k)
      plot(x(:,k) ,zdum(:,i))
      title(['Input ',num2str(k),'-> Output ',num2str(i)]);
    end
  end
end 

thl=[thl(:);shift';scale'];
return;












