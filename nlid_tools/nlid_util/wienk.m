function [k0,k1,k2]=wienk(x,y,numlags,numsides)
% WIENK - compute the 0,1,AND 2 Wiener kernels of a system 
%	[k0,k1,k2] OF A GIVEN SYSTEM WITH THE SPECIFIED INPUT x AND 
%	OUTPUT y.  THE REPEATED TOEPLITZ INVERSION METHOD IS USED.
%
%	Usage :[k0,k1,k2]=wienk(x,y,numpts,numlags,numsides)
%
%	x	: system input
%	y	: system output
%
%	numlags : the number of lags to calculate for both the
%		  first and second order kernels.
%	numsides: determine a causal (1 side) or noncausal 
%		  (2 sides) response.  The noncausal response is
%		  computed by shifting the input forward in time
%		  by numlags/2 points, hence, for a two sided
%		  identification, the following condition must
%		  hold	numpts < size(x) - numlags/2
%
% $Revision: 1.3 $
% EJP Jan 1991
% REK June 1994


% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 



numpts=length(y);
if numsides == 2
     halflen = round(numlags/2);	
     x = x(halflen:numpts + halflen);
end
[rx,cx]=size(x);
[ry,cy]=size(y);  
if cx > 1
     x=x';
end
if cy > 1
     y=y';
end   
x = x - mean(x);    
k0=mean(y);
y = y - k0;
Z=nldat(cat(2,x,y));
I=irf(Z,'nLags',numlags,'nSides',1);
k1=double(I);

y=y-filter(k1,1,x);

x = nldat(x);
phxx = cor(x,'kernOrder',1,'nLags',numlags,'biasMode','biased',...
    'corType','covar','nSides',1);

Z=nldat(cat(2,x,y));
ph2 = cor(phxx,Z,'kernOrder',2);


phxx = double(phxx);
ph2 = double(ph2);

for i=1:numlags
	g(i,:)=0.5*toep(phxx,ph2(i,:))';
end
for i=1:numlags
	k2(:,i)=toep(phxx,g(:,i));
end

