function y = squeeze(x);
% squeeze for nldat objects
y=x;
D=squeeze(double(x));
set(y,'Data',D);
y.Comment=[x.Comment '; squeeze'];


% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
