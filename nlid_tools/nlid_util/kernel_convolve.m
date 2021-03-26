function kernels = kernel_convolve(kerns,subsys)
% KERNEL_CONVOLE - convolve the kernels in kerns with the irf or kernel object in subsys
% kerns is assumed to be a Q element cell array containing kernels of order
% 0 through Q-1
% subsys is assmed to be an irf or volterra kernel.
% $Revision: 1.3 $

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 



ss = get(subsys,'dataSet');
q = get(subsys,'kernOrder');
Ts = get(subsys,'domainIncr');

switch q
  case 0
    % convolve with a zero order kernel, so previous elements are not
    % relevant. 
    kold = kerns;
    Q = length(kold);
    for i = 1:Q
      kold{i} = 0*kold{i};
    end
    kold{1} = subsys;
%    kernels = kerns;
%    set(kernels,'elements',kold);
    kernels = k_old;
  case 1
    k_old = kerns;
    Q = length(k_old);
    for i = 1:Q
      tensor = get(k_old{i},'dataSet');
      kernel = ss(1) * tensor;
      for j = 1:length(ss)-1
        tensor = kernel_delay(tensor);
        kernel = kernel + tensor * ss(j+1);
      end
      temp = k_old{i};
      set(temp,'dataSet',Ts*kernel);
      k_old{i} = temp;
    end
    kernels = k_old;
%    keyboard
%    set(kernels,'elements',k_old);
  otherwise 
    error('high order Volterra kernel in middle of a path');
end
return

function kernel2 = kernel_delay(kernel);
% takes a kernel of arbitrary order (>0), and delays it by one lag

kernel2 = zeros(size(kernel));
hlens = size(kernel);
hlen = hlens(1);
order = length(hlens);

if hlen == 1
  % we have a zero order kernel, so there is no convolution, hence no delay.
  kernel2 = kernel;
elseif min(hlens) == 1  
  % first order kernel, so delay it
  kernel2(2:hlen) = kernel(1:hlen-1);
else
  command = 'kernel2(';
  for i = 1:order-1
    command = [command,'2:hlen,'];
  end
  command = [command,'2:hlen) = kernel('];
  for i = 1:order-1
    command = [command,'1:hlen-1,'];
  end
  command = [command,'1:hlen-1);'];
  eval(command);
end
return
