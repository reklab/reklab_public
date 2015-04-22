function y = nlsim ( ws, u )
% Simulate response of a Wiener series model  to an input signal
% input options not fully defined as yet

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% subtract the input mean, since Wiener series assume a 
% zero mean input.
u = u - mean(u);

% convert the Wiener series into a Volterra series
vs = vseries(ws);

% filter using the Volterra series.
y = nlsim(vs,u);

% ... wseries/nlsim
