clear all
n=4;l=3;m=2;N=100;Ts=0.01;
[A,B,C,D]=rmodel(n,l,m);
u=randn(N,m);
x0=rand(n,1);
K=rand(n,l);
y=tlsim(A,B,C,D,u,Ts,x0);

partype='on'
[theta,params,T]=css2th(A,C);
[Ar,Cr]=cth2ss(theta,params);
T*Ar*inv(T)-A
Cr*inv(T)-C
partype='on'
[theta,params,T]=css2th(A,B,C,D,x0,K,partype);
%reconstruct
[Ar,Br,Cr,Dr,x0r]=cth2ss(theta,params);
yr=tlsim(Ar,Br,Cr,Dr,u,Ts,x0r);
V=vaf(y,yr)

partype='tr'
[theta,params]=css2th(A,B,C,D,x0,K,partype);
[Ar,Br,Cr,Dr,x0r]=cth2ss(theta,params);
yr=tlsim(Ar,Br,Cr,Dr,u,Ts,x0r);
V=vaf(y,yr)
