function coeffs = scale_hermite(mh,sigma);
% produces the normalized hermite coefficients, given a hermite polynomial,
% and the desired input standard deviation.

mhs = mh;
coeffs = get(mh,'Coeffs');
q = length(coeffs);
order = q - 1;
for i = 1:order
  coeffs(i+1) = coeffs(i+1)/(sigma^i);
end

% Copyright 1998-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 
