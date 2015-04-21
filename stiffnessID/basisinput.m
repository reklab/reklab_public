% change the input to satisfy the requirement.
function [U ] = basisinput(u,n)

% this function is to generate the Chebshev polynomial terms 
% of the input that is used to identify Hammerstein system.
% u is the input signal
% y is the constructed inputs
% n is the order of the Chebshev polynomical



if nargin ==0
    help basis_input
end

if size(u,2)>size(u,1)
    u = u';
end

nrow = size(u,1); % length of the data
l = size(u,2); % number of the input

U = zeros(nrow,n);

% construce chebchev polynomial
 % the fist two column

shift=(max(u)+min(u))/2;
scale=(max(u)-min(u))/2;
%scale=(max(u)-min(u))/2;
u=(u-ones(nrow,l)*shift)./(ones(nrow,l)*scale);

U = [ones(size(u)) u];
% u =u;
%U = 0.02*U;
%for j=3:n+1
%    U(:,(j-1)*l+1:j*l)=2*u.*U(:,(j-2)*l+1:(j-1)*l)-U(:,(j-3)*l+1:(j-2)*l);
%end

for j=3:n+1
    U(:,(j-1)*l+1:j*l)=2*u.*U(:,(j-2)*l+1:(j-1)*l)-U(:,(j-3)*l+1:(j-2)*l);
end
U=U(:,2:(n+1)*l);
%U=U';
shift_scale = [shift; scale];