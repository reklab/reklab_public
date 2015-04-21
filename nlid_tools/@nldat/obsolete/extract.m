function y = extract (x, n1, n2);
% Extract samples from a channel dropping first n1 and last n2 points

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if nargin==2,
    n2=n1;
end

[nsamp, nchan, nreal]=size (x);
i1 = n1+1;
i2 = nsamp - n2;
y=x;
dx=double(x);
y.dataSet = dx(i1:i2,:,:);
set (y,'comment','Extract');
return
