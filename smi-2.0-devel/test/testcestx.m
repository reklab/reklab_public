n=4;l=3;m=2;N=100;Ts=0.01;
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
x0=randn(n,1);
y=tlsim(A,B,C,D,u,Ts,x0);
x0e=cestx(u,y,Ts,A,B,C,D)
ye=tlsim(A,B,C,D,u,Ts,x0e);
vaf(y,ye)

