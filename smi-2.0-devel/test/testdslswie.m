clear all
n=4;l=3;m=2;N=100;nn=7;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);
zsim = (y<=.5)*.5+((y>.5)&(y<1.5)).*y+(y>=1.5)*1.5;
z = zsim+.05*randn(size(zsim))*diag(std(zsim));

[thl,ze1]=chebest(y,z,nn);
V1=vaf(z,ze1)

[Ae,Be,Ce,De,x0,thle] = dslswie(u,z,A,B,C,D,[],nn,[],[],1);
ye=dlsim(Ae,Be,Ce,De,u);
ze2=chebsim(ye,thle);
V2=vaf(z,ze2)

if min(V2)>90 ,disp('test succeeded'),else disp('test failed');V,end

[Ae,Be,Ce,De,x0,thle] = dslswie(u,z,A,B,C,[],[],nn,[1,0,0],[],0);
[Ae,Be,Ce,De,x0,thle] = dslswie(u,z,A,[],C,D,[],nn,[0,1,0],[],0);
[Ae,Be,Ce,De,x0,thle] = dslswie(u,z,A,[],C,[],[],nn,[0,0,0],[],0);

disp('test succeeded')
