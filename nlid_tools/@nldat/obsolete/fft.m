function z = fft(x);
% fft function for nldat variables;
% z = fft(x);
%
% matlab's fft is applied to each realization of all channels in data set
% $Revision: 1.1 $


% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

[nsamp, nchan, nreal]=size(x);
z=x;
z.domainIncr= 1/(nsamp*x.domainIncr);
set(z,'comment','FFT');
z.domainStart=0;
z.domainName='Frequency Hz'; 
for ichan=1:nchan,
    for ireal=1:nreal,
        z.dataSet(1:end,ichan,ireal)=fft(x.dataSet(1:end,ichan,ireal));
    end
end
