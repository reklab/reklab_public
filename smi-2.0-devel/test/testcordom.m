n=4;l=3;m=2;N=100;t=0.01;a=20;i=7;
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
y=tlsim(A,B,C,D,u,t);


[S,R]=cordom(u,y,t,a,i);
[S,R2]=cordom(u,y,t,a,i,R);
S
semilogy(S,'x')

