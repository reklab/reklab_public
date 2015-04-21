function [kernel] = lnlkernel ( g,h,Q,kern_len)
%  computes the Q'th order Volterra kernel of an LNL system
%
%  Syntax: kernel = lnlkernel ( g,h,Q,kern_len)
%
%  creates the order kernel for an LNL  cascade with linear elements
%  g followed a squarer followed by h
%
%  kern_len is the maximum size of the kernel
%

% Copyright 1992-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

g = g(:);
h = h (:);
hlen = length(h);
glen = length(g);
if nargin == 2
  kern_len = hlen + glen;
 end

if hlen < kern_len
  h = [h;zeros(kern_len-hlen,1)];
 else
  h = h(1:kern_len);
 end

if glen < kern_len
  g = [g;zeros(kern_len-glen,1)];
 else
  g = g(1:kern_len);
 end

  
tensor = wskern(g,Q);
kernel = h(1) * tensor;
for i = 1:kern_len -1
  tensor = kernel_delay(tensor);
  kernel = kernel + tensor * h(i+1);
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function kernel2 = kernel_delay(kernel);
% takes a kernel of arbitrary order (>0), and delays it by one lag

kernel2 = zeros(size(kernel));
hlens = size(kernel);
hlen = hlens(1);
order = length(hlens);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





