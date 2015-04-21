function z = angle(x);
% angle function for nldat variables;
%

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


z=x;
z.Data=angle(x.Data);
z.ChanUnits='radians';
z.Comment=[x.Comment  ';angle'];
