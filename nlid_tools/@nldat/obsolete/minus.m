function z = minus (x,y);
% minus function for nldat variables;

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

x=nldat(x);
y=nldat(y);
sx=size(x);
xns =sx(1);
xchan=sx(2);
xreal=sx(3);
sy= size(y);
yns=sy(1);
ychan=sy(2);
yreal=sy(3);
s1= [ 1 1 1];
if sx==sy 
   z=x;
   z.Data= x.Data - y.Data;
elseif sx == s1 
   z = y;
   z.Data= x.Data - y.Data;
elseif sy == s1,
   z = x;
   z.Data= x.Data - y.Data;
elseif  sx == [1 ychan 1 ]
   z=y;
   for i=1:xchan,
      z.Data(:,i,:) = x.Data(:,i,:) - y.Data(:,i,:);
   end
elseif sy == [1 xchan 1 ];
   z=x;
   for i=1:xchan,
      z.Data(:,i,:) = x.Data(:,i,:) - y.Data(:,i,:);
   end
   % Separate value for each realization
elseif sy == [1 xchan xreal ];
   z=x;
   for i=1:xchan,
      for j=1:xreal,
         z.Data(:,i,j) = x.Data(:,i,j) - y.Data(:,i,j);
      end  
   end
% subtracting same value from all realizations 
   elseif sy == [xns xchan 1 ];
   z=x;
   for i=1:xchan,
      for j=1:xreal,
         z.Data(:,i,j) = x.Data(:,i,j) - y.Data(:,i,1);
      end  
   end
else
   error ('Dimension missmatch');
end
