function WK = wkkern (z, varargin)
% Wiener kernel model
%  Parent: kern

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

WK=mkwkern;
if nargin==0;
   return
elseif nargin==1,
   WK=nlmkobj(WK,z);
else
   args=varargin;
   WK=nlmkobj(WK,z,args);
end

function w =mkwkern;
% Structure of wkern
w.ObjName='Wiener kernel';
k=kern;
w=class(w,'wkern',k);
w=pdefault(w);
return

% .../@wkern/vern
