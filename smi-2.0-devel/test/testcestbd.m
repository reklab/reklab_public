clear all
n=4;l=3;m=2;N=100;Ts=0.01;
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
y=tlsim(A,B,C,D,u,Ts);
[Be,De]=cestbd(u,y,Ts,A,C);
ye=tlsim(A,Be,C,De,u,Ts);
vaf(y,ye)

x0=randn(n,1);
y=tlsim(A,B,C,D,u,Ts,x0);
[Be,De,x0e]=cestbd(u,y,Ts,A,C,[1 1 1]);
ye=tlsim(A,Be,C,De,u,Ts,x0e);
vaf(y,ye)

[Be,De,R]=cestbd(u,y,Ts,A,C);
[Be,De,R,Phi]=cestbd(u,y,Ts,A,C,[0 1 0]);



[Be,De,xe,Ke,R]=cestbd(u,y,Ts,A,C,[1 1 1 1]);
