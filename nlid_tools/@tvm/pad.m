function Nout = pad (Nin, Ntop, Nbottom); 
% Pad top and bottom of tvm model

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

if nargin==1,
    Nbottom=0;
end


hd=get(Nin,'data');
hfirst = hd{1};
hlast= hd{length(hd)};
i=1;
ds=get(Nin,'DomainStart');
dt=get(Nin,'DomainIncr');
j=1;
hnew={};
for i=1:Ntop,
    hnew{j,1}=hfirst;
    ds=ds-dt;
    j=j+1;
end
for i=1:length(hd),
    hnew{j,1}=hd{i};
      j=j+1;
end
for i=1:Nbottom,
     hnew{j,1}=hlast;
    j=j+1;
end
Nout=Nin;
set(Nout,'data',hnew,'DomainStart',chop(ds,4));
