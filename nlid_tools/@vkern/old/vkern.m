function VK = wvkern (z, varargin)
% Volterra kernel object
%  Parent: kern

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

VK=mkvkern;
if nargin==0;
   return
elseif nargin==1,
  VK=nlmkobj(VK,z);
   
else
   args=varargin;
   VK=nlmkobj(VK,z,args);
end


function v =mkvkern;
v.ObjName='Volterra Kernel';
K=kern;
v=class(v,'vkern',K);
v=pdefault(v);
return

% @vkern/vkern
