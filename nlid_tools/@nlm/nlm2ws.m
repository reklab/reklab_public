function ws = nlm2ws(nlmodel, wsin)
% generate Wiener kernels for a NLM model

% conver the NLM into a Volterra series,
vs = vseries(nlmodel);

% Conver the volterra series into a Wiener series.
ws = wseries(vs,wsin);

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
