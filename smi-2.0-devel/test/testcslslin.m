clear all
n=4;l=3;m=2;N=100;Ts=0.01
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
y=tlsim(A,B,C,D,u,Ts);


[Ae,Be,Ce,De]=cslslin(u,y,Ts,A,C);
ye=tlsim(Ae,Be,Ce,De,u,Ts);
V=vaf(y,ye)


N=1000;
u=rand(N,m);
Q=eye(n);
R=eye(l);
G=eye(n);
v=rand(N,l);
w=rand(N,n);
y = tlsim(A,[B G],C,[D zeros(l,n)],[u, w],Ts);

%theoretical kalman filter
L = dlqe(A,G,C,Q,R);
K = A*L;
yk = tlsim(A-K*C,[B-K*D,K],C,[D zeros(l)],[u y],Ts);
V1=vaf(y,yk)

[Ae,Be,Ce,De,x0e,Ke] = cslslin(u,y,Ts,A,C,K,[1 1 0 1],[],1);
yke = tlsim(Ae-Ke*Ce,[Be-Ke*De,Ke],Ce,[De zeros(l)],[u y],Ts);
V2=vaf(y,yke)








