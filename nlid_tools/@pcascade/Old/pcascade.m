function PC = pcascade (z, varargin)
% Parallel Cascade Model
%  Parent: nlm

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

PC=mkpc;
if nargin==0;
   return
elseif nargin==1,
   PC=nlmkobj(PC,z);
else
   args=varargin;
   PC=nlmkobj(PC,z,args);
end

    
function pc = mkpc
S=nlm;
i=irf;
t=polynom;
set(t,'polyType','tcheb');
elements = {i t; i t};
set (S,'Elements',elements);
pc.Type='Parallel Cascade';
pc=class(pc,'pcascade',S);
set(pc,'Method','eig');
pc=pdefault(pc);



return


