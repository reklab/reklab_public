function wkernel = lnl2wiener(lnl,hlen,order,sigma_u);
% compute first and second-order Wiener kernels of an LNL system
%
% syntax:  wkernel = lnl2wiener(lnl,hlen,order,sigma_u);

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


% still to do -- non-zero mean inputs 


% unpack the lnl object.

elements = get(lnl,'elements');
Ts = get(elements{1},'domainincr');
g = get(elements{1},'data');
m = elements{2};
max_order = get(m,'Order');
h = get(elements{3},'data');

% compute the Std of the x, the output of g, the first linear subsystem. 
sigma_x = Ts * sigma_u * sqrt(sum(g.^2));


% convert the nonlinearity into a hermite polynomial, normalized for an
% input with Std sigma_x
polytype = get(m,'type');
if strcmp(polytype,'hermite')
  m = polynom(m,'type','tcheb');
end
set(m,'Std',sigma_x,'Mean',0);
mh = polynom(m,'type','hermite');
coeffs = get(mh,'Coef');


wkernel = wkern;
set(wkernel,'order',order,'nlags',hlen,'domainincr',Ts);


if order > max_order
  set(wkernel,'data',0,'nlags',1);
elseif order < 0
  error('order must be positive');
else
  switch order
    case 0
      % zero-order Wiener kernel is the DC gain of h multiplied by the
      %  zero-order hermite polynomial coefficient of the nonlinearity
      set(wkernel,'nlags',1);
      kernel = Ts*sum(h)*coeffs(1);
      set(wkernel,'data',kernel);
    case 1
      % first-order Volterra kernel is g*h times the first-order normalized
      % hermite polynomial coefficient.
      kernel = (coeffs(2)/sigma_x)*Ts*conv(g,h);
      if length(kernel) > hlen
         kernel = kernel(1:hlen);
      else
        pad = zeros(hlen - length(kernel),1);
        kernel = [kernel; pad];
      end
      set(wkernel,'data',kernel);
    otherwise
      % order 'order' Volterra kernel is computed using lnlkernel (below)
      kernel = (coeffs(order+1)/(sigma_x^order)) * Ts *...
	  lnlkernel(g,h,order,hlen);
      set(wkernel,'data',kernel);
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
%  DTW Jan 1992
%  DTW Jan 1995  modified handling of IRF lengths and kernel size.
%
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

  
tensor = wsirf2kern(g,Q);
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




