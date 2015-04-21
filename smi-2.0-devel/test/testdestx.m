n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
x0=randn(n,1);
y=dlsim(A,B,C,D,u,x0);
x0e=destx(u,y,A,B,C,D);
ye=dlsim(A,B,C,D,u,x0e);
vaf(y,ye)

