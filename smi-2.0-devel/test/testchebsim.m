n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);
zsim = (y<=.5)*.5+((y>.5)&(y<1.5)).*y+(y>=1.5)*1.5;
z = zsim+.05*randn(size(zsim))*diag(std(zsim));

[thl]=chebest(y,z,7);
ze=chebsim(y,thl);
vaf(z,ze)

[thl]=chebest(y,z,8);
ze=chebsim(y,thl);
vaf(z,ze)
