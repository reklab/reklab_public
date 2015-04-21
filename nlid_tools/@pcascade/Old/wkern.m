function ws = wseries (pc);
% Generate serioes from parallel cascade description
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

num_paths=size(pc);
hlen=get(pc{1,1},'NLags');
k0 = 0;
k1 = zeros(hlen,1);
k2 = zeros(hlen,hlen);
for i=1:num_paths,
   h =double (pc{i,1});
   p=set(pc{i,2},'type','power');
   mc=get(p,'Coef');
   k0 = k0 + mc(1);
   if length (mc) > 1,
      k1 = k1 + mc(2)*h;
      if length(mc) >2,
         k2 = k2 + mc(3)*h*h';
      end
   end
end
k0=wkern(k0);
k1=wkern(k1);
k2=wkern(k2);
k = { k0; k1; k2 };
wk=wseries;
set(ws,'elements',k);





  
