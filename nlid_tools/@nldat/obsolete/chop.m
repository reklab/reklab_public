function z = chop (x, chopLen);
% nldat/chop - chop data from begining and end of a record 
% z = chop (x, chopLen);
% x - input data
% chop len - length tyo chop off
% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


istart = round( (chopLen)/x.domainIncr) + 1 ;
if istart>length(x.dataSet/2),
    error ('chopLen too large');
end
iend= length(x.dataSet) -istart + 1;
x.dataSet=x.dataSet(istart:iend,:,:);
x.domainStart=x.domainStart+chopLen;
Comment = x.comment;
set (x,'comment',[ Comment '; Chop']);
z=x;
