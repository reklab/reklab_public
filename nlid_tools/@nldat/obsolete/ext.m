function z = ext (x, start, delt);
% extract data based on domain values
% x - input data
% start - starting domain vlaue 
% delt  - length of extract. (default is to end of data set);
% z = x(start:start_delt, :,:);

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


istart = (start-x.domainStart)/x.domainIncr; 
if nargin == 2,
  iend=length(d.dataSet);
else
  iend=min(istart + delt/x.domainIncr, length(x.dataSet));
end
x.dataSet=x.dataSet(istart:iend,:,:);
x.domainStart=start;
comment = get (x,'comment');
set (x,'comment',[ comment '; ext']);
z=x;
% nldat/ext
