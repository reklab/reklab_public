function y = mean (x, DIM);
y=nldat(x);
if nargin==1,
  DIM=1;
end
set(y,'Comment',['Mean(' inputname(1) ',' int2str(DIM) ')' ]);
y.Data=mean(x.Data,DIM);
%     if DIM==1,
%        [nx,ny,nz]=size(y.Data);
%        z=reshape(y.Data,ny,nz);
%        y.Data=z';
%        y.DomainIncr=1;
%        y.DomainStart=1;
%        y.DomainName='Realization';
%        for i=1:ny,
%           y.ChanNames{i}=['mean(' y.ChanNames{i} ')'];
%        end    
%     end

return
% nldat/mean

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
