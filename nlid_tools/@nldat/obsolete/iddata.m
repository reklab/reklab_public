function xid = iddata (N)
% converts nldata type to iddata type
% only for single-input singple output

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

xin=N.Data(:,1);
xout=N.Data(:,2);
xid=iddata(N.Data(:,2),N.Data(:,1));
set(xid,'Ts',N.DomainIncr');
return
