function X = flb2nld ( f, N ) 
% FLB2NLD - Read a nld object from flb
% $Revision: 1.5 $
% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if nargin == 1,
   x=flb2mat(f);
   return
end
% Read a nldat type from a flb file 
x= flb2mat(f,'read_Case',N);
X=nldat;
set(X,'dataSet', x.Data,'comment',x.comment,'domainIncr', x.domainIncr,'domainName',x.domainName,'chanNames',x.chanName);
end