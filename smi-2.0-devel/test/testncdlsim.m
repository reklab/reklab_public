n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);

while (max(abs(eig(A)))<1)
  A=A*1.1;
end
y=ncdlsim(A,B,C,D,u);

n=4;l=1;m=1;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);

while (max(abs(eig(A)))<1)
  A=A*1.1;
end
y1=ncdlsim(A,B,C,D,u);
[num,den]=ss2tf(A,B,C,D);
y2=ncdlsim(num,den,u);
vaf(y1,y2)


[A,B,C,D]=drmodel(n,l,m);
[Ae,Be,Ce]=drmodel(n,l,m);
while (max(abs(eig(Ae)))>1)
  Ae=Ae*0.9;
end
y1=ncdlsim(A,Ae,B,Be,C,Ce,D,u);
