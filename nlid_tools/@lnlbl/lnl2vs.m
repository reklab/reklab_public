function vs=lnl2vs (lnl, vsin)
% generate Volterra kernels for a LNL cascade model

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

subsys = get(lnl,'elements');

h = double(subsys{1});
hlen = length(h);
Ts = get(subsys{1},'domainincr');

p = subsys{2};
p = nlident(p,'type','power');
mc=get(p,'Coef');
Q = length(mc);
kernels = cell(Q,1);

g = double(subsys{3});
glen = length(g);


% If the memory length of the Volterra kernels was not specified, use the
% total memory length of the LNL cascade.  Otherwise, use the memory length
% specified by vsin.  
NlagsVs = get(vsin,'nlags')
if isnan(NlagsVs)
  % use total LNL memory length -- pad h with zeros.
  hpad = [h;zeros(glen,1)];
elseif NlagsVs > hlen
  % pad h with zeros out to the kernel length
  hpad = zeros(NlagsVs,1);
  hpad(1:hlen) = h;
else
  % truncate h to the kernel length.
  hpad = h(1:NlagsVs);
end

  
% generate the kernels of the first two elements.

irf1 = subsys{1};
set(irf1,'data',hpad);

ln = lnbl;
set(ln,'elements',{irf1, p});
vs = ln2vs(ln,vsin);
kernels = get(vs,'elements');

% convolve the LN kernels with the second linear element.
kernels = kernel_convolve(kernels,subsys{3});

set(vs,'elements',kernels,'comment','Transformed LNL Cascade');


