function plot (vs)
% Plot a volterra series model
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

kernels =get(vs,'Elements');
Qmax = length(kernels)-1;

switch Qmax
  case 0
    subplot(111);
    plot(kernels{1});
  case 1
    subplot(211)
    plot(kernels{1});
    subplot(212)
    plot(kernels{2});
  case 2
    subplot(221)
    plot(kernels{1});
    subplot(222)
    plot(kernels{2});
    subplot(212)
    plot(kernels{3});
  otherwise
    subplot(221)
    plot(kernels{1});
    subplot(222)
    plot(kernels{2});
    subplot(223)
    plot(kernels{3});
    subplot(224)
    plot(kernels{4});
end


streamer (get(vs,'Comment'));


