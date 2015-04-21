function z = times(x,y);
% array multiply function for nldat variables;

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

x=nldat(x);
y=nldat(y);
sx=size(x);
xchan=sx(2);
sy= size(y);
ychan=sy(2);
s1= [ 1 1 1];
if sx==sy 
   z=x;
   z.Data= x.Data .* y.Data;
elseif sx == s1 
   z = y;
   z.Data= x.Data .* y.Data;
elseif sy == s1,
   z = x;
   z.Data= x.Data .* y.Data;
elseif  sx == [1 ychan 1 ]
   z=y;
   for i=1:xchan,
      z.Data(:,i,:) = x.Data(:,i,:) .* y.Data(:,i,:);
   end
elseif sy == [1 xchan 1 ];
   z=x;
   for i=1:xchan,
      z.Data(:,i,:) = x.Data(:,i,:) .* y.Data(:,i,:);
   end
else
   error ('Dimension missmatch');
end
