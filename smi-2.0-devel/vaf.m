function v = vaf(y,ye)
% vaf        Compute the percentage Variance Accounted For (VAF)
%            between two signals. the VAF is calculated as:
%                         variance(y-ye)
%               v= ( 1 -  ----------------  ) * 100%
%		            variance(y)
% 
%            The VAF of two signals that are the same is
%            100%. If they differ, the  VAF will be lower. 
%            When y and ye have multiple columns, the VAF 
%            is calculated for every column in y and ye. 
%            The VAF is often used to verify the
%            correctness of a model, by comparing the real
%            output with the estimated output of the model.
% 
% Syntax:
%            v = vaf(y,ye)
% Input:     
%  y         Reference signal, often the real output.
%  ye        Second signal, often the estimated output of a model.
%
% Output:
%   v        VAF, computed for the two signals

% Bert Haverkamp, April 1996
% copyright (c) 1996 B.R.J. Haverkamp

if nargin==0
  help vaf
  return
end
if nargin<2
  error('Not enough input variables')
end
if size(y,2)>size(y,1)
  y=y';
end
if size(ye,2)>size(ye,1)
  ye=ye';
end

N=size(y,1);

if ~(size(ye,1)==N)
  error('Both signals should have same lenght')
end

v=max(diag(100*(eye(size(y,2))-cov(y-ye)./cov(y))),-1000);





