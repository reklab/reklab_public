function [X,r] = pinv(A,tol)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS.
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%
%   See also RANK.

%   Copyright (c) 1984-98 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2008-10-21 14:28:03 $

%
% This function now also returns the number of singular vectors used to compute the pseudoinverse.
% Mireille 1998.09.27
%

[U,S,V] = svd(A,0);
[m,n] = size(A);
if m > 1, s = diag(S);
   elseif m == 1, s = S(1);
   else s = 0;
end
if nargin < 2
   tol = max(m,n) * max(s) * eps;
end
r = sum(s > tol);
if (r == 0)
   X = zeros(size(A'));
else
   s = diag(ones(r,1)./s(1:r));
   X = V(:,1:r)*s*U(:,1:r)';
end
