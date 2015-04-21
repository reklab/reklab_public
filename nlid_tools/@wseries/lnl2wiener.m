function wkernel = lnl2wiener(lnl,hlen,order,sigma_u);
% compute first and second-order Wiener kernels of an LNL system
%
% syntax:  wkernel = lnl2wiener(lnl,hlen,order,sigma_u);

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% still to do -- non-zero mean inputs 


% unpack the lnl object.

elements = get(lnl,'elements');
Ts = get(lnl,'Ts');
g = get(elements{1},'Kernel');
m = elements{2};
max_order = get(m,'Order');
h = get(elements{3},'Kernel');

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
coeffs = get(mh,'Coeffs');


wkernel = wkern;
set(wkernel,'order',order,'nlags',hlen,'Ts',Ts);


if order > max_order
  set(wkernel,'Kernel',0,'nlags',1);
elseif order < 0
  error('order must be positive');
else
  switch order
    case 0
      % zero-order Wiener kernel is the DC gain of h multiplied by the
      %  zero-order hermite polynomial coefficient of the nonlinearity
      set(wkernel,'nlags',1);
      kernel = Ts*sum(h)*coeffs(1);
      set(wkernel,'Kernel',kernel);
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
      set(wkernel,'Kernel',kernel);
    otherwise
      % order 'order' Volterra kernel is computed using lnlkernel (below)
      kernel = (coeffs(order+1)/(sigma_x^order)) * Ts *...
	  lnlkernel(g,h,order,hlen);
      set(wkernel,'Kernel',kernel);
  end      
end




