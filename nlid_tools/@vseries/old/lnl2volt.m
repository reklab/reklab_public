function VKernel = lnl2volt(lnl,hlen,order);
%  computes Volterra kernels of a lnl cascade.
%
%  syntax  VKernel = lnl2volt(lnl,hlen,order);

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

elements = get(lnl,'elements');
Ts = get(lnl,'Ts');
g = get(elements{1},'Kernel');
m = elements{2};
mc = polynom(m,'type','power');
coeffs = get(mc,'Coeffs');
max_order = get(m,'Order');
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



    
