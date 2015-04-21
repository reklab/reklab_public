function y = seg (x,N);
% segment a nldat object into a series of relizations N points long
[nx,ny,nz]=size(x);
nseg = nx/N;
y=x;
j=1;
y.Data=[];
for iz=1:nz,
   
   for i=1:N:nx
      y.Data(1:N,:,j)=x.Data(i:(i+N-1),:,iz);
      j=j+1;
   end
end

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
