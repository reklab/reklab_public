clear all
n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
x0=rand(n,1);
K=rand(n,l);
u=randn(N,m);
y=dlsim(A,B,C,D,u);

partype='tr';
[theta,params]=dss2th(A,B,C,D,x0,K,partype);
%reconstruct
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
yr=dlsim(Ar,Br,Cr,Dr,u);
V=vaf(y,yr);
if V==100,disp('test succeeded'),else disp('test failed');V,end

partype='on';
[theta,params]=dss2th(A,B,C,D,x0,K,partype);
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
yr=dlsim(Ar,Br,Cr,Dr,u);
V=vaf(y,yr);
if V==100,disp('test succeeded'),else disp('test failed');V,end

[theta,params]=dss2th(A,[],C,D,x0,K,partype);
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
[theta,params]=dss2th(A,B,C,[],x0,K,partype);
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
[theta,params]=dss2th(A,B,C,D,[],K,partype);
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
[theta,params]=dss2th(A,[],C,[],[],K,partype);
[Ar,Br,Cr,Dr]=dth2ss(theta,params);
disp('Test succeeded')




