clear all
n=4;l=3;m=2;N=100;Ts=0.01;
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
x0=rand(n,1);
K=rand(n,l);
y=tlsim(A,B,C,D,u,Ts,x0);



partype='on';
theta=css2th(A,C)
theta=css2th(A,C,partype)
theta=css2th(A,B,C,D)
theta=css2th(A,B,C,D,partype)
theta=css2th(A,B,C,D,x0,partype)
theta=css2th(A,B,C,D,x0,K,partype)


partype='tr';
theta=css2th(A,C)
theta=css2th(A,C,partype)
theta=css2th(A,B,C,D)
theta=css2th(A,B,C,D,partype)
theta=css2th(A,B,C,D,x0,partype)
theta=css2th(A,B,C,D,x0,K,partype)
