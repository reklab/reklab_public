function wskernel = pcas2wiener(pc_system,order,sigma_u);
% compute first and second-order Wiener kernels of an LNL system
%
% syntax:  wskernel = lnl2wiener(pc_model,order,sigma_u);

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% still to do -- non-zero mean inputs 


% 1. decompose the Parallel Cascade model, convert the Polynomials into
%    Normalized Hermite form.

num_paths = get(pc_system,'NPaths');
if order >0
   hlen = get(pc_system,'NLags');
else
   hlen=0;
end

Ts =  get(pc_system,'Ts');
max_order =  get(pc_system,'POrderMax');
elements = get(pc_system,'elements');


pcH = zeros(hlen,num_paths);
pcM = zeros(max_order+1,num_paths);

for i = 1:num_paths
   
   h = get(elements{i,1},'Kernel');
   lh=length(h);
   pcH(1:lh,i) = h;
   nh_poly = NormHermite(elements{i,1}, sigma_u, elements{i,2});
   p = get(nh_poly,'Coeffs');
   lp=length(p);
   pcM(1:lp,i) = p;
end

switch order
  case 0
    kernel = sum(pcM(1,:));
  case 1
    kernel = pcH * pcM(2,:)';
  otherwise
    kernel =  pcM(order+1,1) * wsirf2kern(pcH(:,1),order);
    for i = 2:num_paths
      kernel =  kernel + pcM(order+1,i) * wsirf2kern(pcH(:,i),order);
    end
end

wskernel = wkern;
set(wskernel,'order',order,'Ts',Ts,'NLags',hlen);
set(wskernel,'kernel',kernel);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nh_poly = NormHermite(h,sigma_u,m)
% converts old_poly into a Normalized Hermite polynomial, 

% compute the Std of the x, the output of h, the linear subsystem. 
Ts = get(h,'Ts');
w = get(h,'Kernel');
sigma_x = Ts * sigma_u * sqrt(sum(w.^2));

% convert the nonlinearity into a hermite polynomial, normalized for an
% input with Std sigma_x
polytype = get(m,'type');
if strcmp(polytype,'hermite')
  m = polynom(m,'type','tcheb');
end
set(m,'Std',sigma_x,'Mean',0);
nh_poly = polynom(m,'type','hermite');
coeffs = get(nh_poly,'Coeffs');
order = get(nh_poly,'Order');
for i = 1:order
  coeffs(i+1) = coeffs(i+1)/(sigma_x^i);
end
set(nh_poly,'Coeffs',coeffs);


