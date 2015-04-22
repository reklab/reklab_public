function y = decimate (d, n, type);
% overloaded decimate fuction for "data" class
% usage  y = decimate (d, n);
% or     y = decimate (d, n, 'fir');

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

% Modified by DL and TSV on 15/9/08 to enable 'fir' filtering

[nsamp,nchan,nreal]=size(d);
y=nldat(d);

if nargin == 3
    if (strcmp(type,'fir'))
        for i=1:nchan,
            for j=1:nreal,
                xin = d.dataSet(:,i,j);
                dout(:,i,j)=decimate(xin,n,'fir');
            end
        end
    else
        error('Only argument is ''fir''')
    end
else
    for i=1:nchan,
        for j=1:nreal,
            xin = d.dataSet(:,i,j);
            dout(:,i,j)=decimate(xin,n);
        end
    end
end
y.dataSet=dout;
y.domainIncr=d.domainIncr*n;
y.domainStart=d.domainStart + (n-1)*d.domainIncr;
y.comment = [ d.comment '; decimate'];

% 24 Jan 2002 Set DomainStart correctly.


