function s=prt (x)
[nrow,ncol,nreal]= size(x);
dv=get(x,'domainvalues');
dt=get(x,'data');
cn=get(x,'channames');
s=sprintf('%8s','-');
for i=1:ncol,
   s=strcat( s,sprintf('%10s',char(cn{i})));
end
s1=[];
for i=1:nreal,
   s1=sprintf ('%8s',char(dv(i)));
   for j= 1:ncol,
      s1=strcat(s1,sprintf ('%10.4f', dt(1,j,i)));
   end
s=strvcat(s,s1);  
end

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
      
