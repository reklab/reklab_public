n=4;l=3;m=2;N=200;i=6;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


[S,R]=dordpi(u,y,2*n);
[Ae,Ce]=destac(R,n);
[Be,De]=destbd(u,y,A,C);
x=ltitr(Ae,Be,u);
[S,R]=dordrs(u,y,x,i);
[S,R2]=dordrs(u,y,x,i,R);
[Ae,Ce]=destac(R,n);
eigA=eig(A),eigAe=eig(Ae)

