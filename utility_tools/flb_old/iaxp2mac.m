function imac = iaxp2mac (iaxp)
%
% Convert a axp integer*4 to a mac integer
%
imac=iaxp(1) + iaxp(2)*(2^8) + iaxp(3)*(2^16) + iaxp(4)*(2^24) ;
return
