function y = nlsim ( vs, u )
% Simulate response of a Volterra series model  to an input signal
% input options not fully defined as yet

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

%Q =get (vs,'order');
kernels = get(vs,'elements');
num_kerns = length(kernels);
ud=double(u);
y=nldat(u);
yd = zeros(size(ud));

for i = 1:num_kerns
  vk = kernels{i};
  yd = yd + double(nlsim(vk,u));
end
set(y,'Data',yd);

% ... vseries/nlsim
