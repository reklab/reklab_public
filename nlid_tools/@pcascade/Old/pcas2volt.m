function vskern = pcas2volt(pc_system,order);
% generates a Volterra kernel of a parallel Wiener cascade model
% 
% syntax:  vskern = pcas2volt(pc_system,order);

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% 1. decompose the Parallel Cascade model, convert the Polynomials into
%    power series form.
% REK modify to handle varying polynominal order
% DTW applied changes to pcascade and vkern object property/parameter names 

num_paths = get(pc_system,'NPaths');
if order > 0
   hlen = get(pc_system,'NLags');
else
   hlen=1;
end

max_order =  get(pc_system,'OrderMax');
elements = get(pc_system,'elements');
Ts = get(elements{1,1},'DomainIncr');

pcH = zeros(hlen,num_paths);
pcM = zeros(max_order+1,num_paths);

for i = 1:num_paths
   e = get(elements{i,1},'Data');
   le=length(e);
   pcH(1:le,i) = e;

  poly = polynom(elements{i,2},'type','power');
  c= get(poly,'Coef');
  lc=length(c);
  pcM(1:lc,i) = c;

end

switch order
  case 0
    vskernel = sum(pcM(1,:));
  case 1
    vskernel = pcH * pcM(2,:)';
  otherwise
    vskernel =  pcM(order+1,1) * wckern(pcH(:,1),order);
    for i = 2:num_paths
      vskernel =  vskernel + pcM(order+1,i) * wckern(pcH(:,i),order);
    end
end

vskern = vkern;
set(vskern,'order',order,'DomainIncr',Ts,'NLags',hlen);
set(vskern,'Data',vskernel);

