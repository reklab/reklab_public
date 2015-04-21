function y = nlsim ( wk, u, uvar )
% Simulate response of wkern to input data set
% input options not fully defined as yet
 
% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

order=get(wk,'order');
Ts=get(wk,'domainincr');
k=get (wk,'data');
ud=double(u(:,1));



if nargin < 3
  warning('orthogonalization variance not specified, using signal');
  uvar = std(ud)^2;
end

% use the Volterra kernel computing routines.
kv = vkern;
set(kv,'domainincr',Ts,'order',order,'data',k);



switch order
case 0
   y=u*0 + k;
 case 1
   y = nlsim(kv,u);
 case 2
   y20 = Ts^2* uvar * sum(diag(k));
   y = nlsim(kv,u)-y20;
 case 3
   k31 = sum(squeeze(sum(k)));
   vk31 = kv;
   set(vk31,'order',1,'data',k31);
   y = nlsim(kv,u)-3*uvar*nlsim(vk31,u);
  otherwise 
    error ('nlsim not defined for kernels of order > 3');
end


return
