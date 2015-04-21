function VKernel = lnl2vs(lnl,hlen,order);
%  computes Volterra kernels of a lnl cascade.
%
%  syntax  VKernel = lnl2volt(lnl,hlen,order);
warning ('iunder development'); 
elements = get(lnl,'elements');
Ts = get(lnl,'Ts');
g = get(elements{1},'Kernel');
m = elements{2};
mc = polynom(m,'polyType','power');
coeffs = get(mc,'polyCoef');
max_order = get(m,'polyOrder');
h = get(elements{3},'Kernel');

VKernel = vkern;
set(VKernel,'order',order,'nlags',hlen,'Ts',Ts);

if order > max_order
  set(VKernel,'Kernel',0,'nlags',1);
elseif order < 0
  error('order must be positive');
else
  switch order
    case 0
      % zero-order Volterra kernel is the DC gain of h multiplied by the
      %  zero-order polynomial coefficient of the nonlinearity
      set(VKernel,'nlags',1);
      kernel = Ts*sum(h)*coeffs(1);
      set(VKernel,'Kernel',kernel);
    case 1
      % first-order Volterra kernel is g*h times the first-order polynomial
      % coefficient.
      kernel = coeffs(2)*Ts*conv(g,h);
      kernel = kernel(1:hlen);
      set(VKernel,'Kernel',kernel);
    otherwise
      % order 'order' Volterra kernel is computed using lnlkernel (below)
      kernel = coeffs(order+1)*Ts*lnlkernel(g,h,order,hlen);
      set(VKernel,'Kernel',kernel);
  end      
end

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









    
