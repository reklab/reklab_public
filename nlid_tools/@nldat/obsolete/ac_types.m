function ac_types

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


function ac_examples
clf;
nx=4;
ny=2;
% {{{  white noise
tmax=2;incr = .001;
x=nlid_sig('normal', 'meanval', 1.,'sd', 3, 'incr', incr, ...
   'tmax',tmax);
stitle='white noise with offset';
n=1;
subplot (nx,ny,n);
plot (x);
axis([-inf inf -10 10]);
fixplot ('A');
acor= xcorr(x,'ctype','corel','maxlag',.01);
acov= xcorr(x,'ctype','cov','maxlag',.01);
acoef= xcorr(x,'ctype','coeff','maxlag',.01);
subplot (nx,ny,n+1);
plot (acor);
subplot (nx,ny,n+2);
plot (acov);
subplot (nx,ny,n+2);
plot (acoef);
fixplot ('B');
ylabel ('cxx');
axis([-inf inf -.1 1.1]);

% }}}
