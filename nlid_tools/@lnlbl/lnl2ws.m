function ws=lnl2ws (lnl, wsin)
% generate Wiener kernels for a LNL cascade model

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


subsys = get(lnl,'elements');

h = double(subsys{1});
hlen = length(h);
Ts = get(subsys{1},'domainincr');

p = subsys{2};
Q = get(p,'order')+1;
kernels = cell(Q,1);

g = double(subsys{3});
glen = length(g);


% Pad the first IRF out to the memory length of the LNL system, 
% and generate the kernels of the first two elements.
hpad = [h;zeros(glen,1)];
irf1 = subsys{1};
set(irf1,'data',hpad);

ln = lnbl;
set(ln,'elements',{irf1, p});
ws = ln2ws(ln,wsin);

%keyboard

kernels = get(ws,'elements');

% convolve the LN kernels with the second linear element.
kernels = kernel_convolve(kernels,subsys{3});

set(ws,'elements',kernels,'comment','Transformed LNL Cascade');
