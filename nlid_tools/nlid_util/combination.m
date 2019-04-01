%% This function calculates the number of k-combinations from n objects where the order of selection does not matter (as opposed to permutation).
function y = combination(n,k)
%++ Instead of computationally in-efficient general formula
% y = factorial(n) / (factorial(n-k)*factorial(k));     

%++ We use the following more efficient implementation: 
a = (n:-1:n-k+1);
y = prod(a) / factorial(k);

end

