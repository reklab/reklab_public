function z = cat (DIM, x, y);
% cat function for nldat objects

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

x=nldat(x);
y=nldat(y);
z=x;
z.Data = cat (DIM, x.Data, y.Data);
switch DIM
case 2
   z.ChanNames=cat(2,x.ChanNames,y.ChanNames);
end
Comment = [ 'cat( ' int2str(DIM) ' ' inputname(2) inputname(3) ];
return
% @nldat/cat
