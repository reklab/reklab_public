function Y = hwrect(X)

Y = X;
x = double(X);
y = (x>0).*x;
Y.Data=y;
Y.Comment = [ X.Comment '; hw rect'];

% nldat/hwrect
% DL 17 Sep  2008

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 