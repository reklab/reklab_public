function VS = vseries (z, varargin)
% Volterra Series Model
%  Parent: nlm

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

VS=mkvseries;
if nargin==0;
   return
elseif nargin==1,
  VS=nlmkobj(VS,z);
else
   args=varargin;
   VS=nlmkobj(VS,z,args);
end


function v =mkvseries;
% Structure of volterra series
N=nlm;
vk=vkern;
v0=vkern(vk,'order',0);
v1=vkern(vk,'order',1);
v2=vkern(vk,'order',2);
el = { v0 ; v1; v2 };
v.ObjType='Volterra Series';
set(N,'Elements',el);
v=class(v,'vseries',N);
set(v,'method','FOA');
v=pdefault(v);
return

% .../@vseries/vseries
