n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);

[Phi,T]=dcalcphi(u,y,A,C,[1,1,1,1,1])