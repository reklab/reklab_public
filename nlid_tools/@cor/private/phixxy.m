function phi=phixxy(x,y,hlen)
%  the second-order cross correlation between  x and y
%
%  syntax: phi = phixxy (x,y,hlen)
%
%  the second-order cross correlation between  x and y
%  output y, out to hlen lags.  
%  the means of both x and y are subtracted prior to estimation
%  a biased estimate of the cross-correlation is returned

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


disp('This function is implemented as a MEX file.');
disp('If MATLAB is giving you this message, the nlid toolbox MEX files');
disp('have not been installed for your machine architecture');

