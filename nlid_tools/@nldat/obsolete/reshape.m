function y = reshape (x, nsamp, nchan, nreal)
% overlaid reshape for nldat objects
if nargin < 4,
    nreal=1;
end
y=nldat(x);
y.ChanNames={};
y.Data = reshape (x.Data, nsamp, nchan, nreal);
y.Size=[ nsamp, nchan nreal];

for i=1:nchan,
    y.ChanNames{i}= 'Reshaped channel';
end

return

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
