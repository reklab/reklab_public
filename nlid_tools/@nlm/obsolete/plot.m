function plot (n)
% Plot a nlm model
%

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

e=n.Elements;
[nout,nin]=size(e);
% SISO Series element
[np,ns]=size(e);
ifig=gcf;
for i=1:np,
   for j=1:ns
      subplot (np,ns,(i-1)*ns+j);
      p=e{i,j};
      plot (p);
      title(n.Comment);
      h=get(gca,'title');
      u=get(h,'units');
      set(h,'units','pixels');
      p=get(h,'position');
      set(h,'position',p+[0 -30 0]);
      set(h,'units',u);
      
   end
end
streamer ((n.Comment));


