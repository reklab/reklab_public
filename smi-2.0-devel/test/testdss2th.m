clear all
n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
x0=rand(n,1);
K=rand(n,l);
u=randn(N,m);
y=dlsim(A,B,C,D,u);

partype='on';
theta=dss2th(A,C)
theta=dss2th(A,C,partype)
theta=dss2th(A,B,C,D)
theta=dss2th(A,B,C,D,partype)
theta=dss2th(A,B,C,D,x0,partype)
theta=dss2th(A,B,C,D,x0,K,partype)
