function z = uminus (x);
% unary function for nldat variables;
% 
z=nldat(x);
z.Data=-z.Data;

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
