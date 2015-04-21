clear all
n=4;l=3;m=2;N=100;Ts=0.01
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
y=tlsim(A,B,C,D,u,Ts);


Ai=A+0.2*randn(size(A));
Ci=C+0.2*randn(size(C));
options=1;
partype='on';
[Ae,Be,Ce,De,hist] = crslslin(u,y,Ts,Ai,Ci,[],[],[],0);
ye1=dlsim(Ai,B,Ci,D,u);
V1=vaf(y,ye1)
ye2=dlsim(Ae,Be,Ce,De,u);
V2=vaf(y,ye2)

