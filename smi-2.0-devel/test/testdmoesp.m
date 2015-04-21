n=4;l=3;m=2;N=100;i=8;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


[Ae,Be,Ce,De]=dmoesp(u,y,i);
ye=dlsim(Ae,Be,Ce,De,u);
V=vaf(y,ye)

[Ae,Be,Ce,De]=dmoesp(u,y,i,n);
ye=dlsim(Ae,Be,Ce,De,u);
V=vaf(y,ye)
[Ae,Be,Ce,De]=dmoesp(u,y,i,n,'pi');
ye=dlsim(Ae,Be,Ce,De,u);
V=vaf(y,ye)