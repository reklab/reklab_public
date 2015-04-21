%make anticausal model 
u=rand(101,1);
y=u(2:101);
u=u(1:100);
[S,R]=dordpo(u,y,4);
[Ace,Aae,Cce,Cae]=ncdestac(R,1);
[Bce,Bae,De,xce,xae]=ncdestbd(u,y,Ace,Aae,Cce,Cae,[1 1 1]);
ye=ncdlsim(Ace,Aae,Bce,Bae,Cce,Cae,De,u,xce,xae);
vaf(y,ye)

n=4;l=3;m=2;N=101;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);
y=u(2:101,:);
u=u(1:100,:);
[S,R]=dordpo(u,y,8);
[Ace,Aae,Cce,Cae]=ncdestac(R,4);
%[Bce,Bae,De]=ncdestbd(u,y,Ace,Aae,Cce,Cae);
[Bce,Bae,De,xce,xae]=ncdestbd(u,y,Ace,Aae,Cce,Cae,[1 1 1]);
ye=ncdlsim(Ace,Aae,Bce,Bae,Cce,Cae,De,u,xce,xae);
vaf(y,ye)

%  test ncdestbdx
nc=3;na=3;n=5;l=2;m=1;N=1000;Np=N/2;
[Ac,Bc,Cc,D]=drmodel(nc,l,m);
[Aa,Ba,Ca]=drmodel(na,l,m);
if max(abs(eig(Ac)))>1
  Ac=Ac*0.95;
end
if max(abs(eig(Aa)))>1
  Aa=Aa*0.95;
end

xc=rand(nc,1);
xa=rand(na,1);
u=randn(N,m);
y=ncdlsim(Ac,Aa,Bc,Ba,Cc,Ca,D,u,xc,xa);
[S,R]=dordpo(u,y,6);
[Ace,Aae,Cce,Cae]=ncdestac(R,n);
[Bce,Bae,De]=ncdestbd(u,y,Ace,Aae,Cce,Cae);
[Bce,Bae,De,xce,xae]=ncdestbd(u,y,Ace,Aae,Cce,Cae,[1 1 1]);
ye=ncdlsim(Ace,Aae,Bce,Bae,Cce,Cae,De,u,xce,xae);
vaf(y,ye)

[Bce,Bae,De,xc0,xa0,R,Phi]=ncdestbd(u,y,Ac,Aa,Cc,Ca,[1,1,0]);
[Bce,Bae,De]=ncdestbd(u,y,Ac,Aa,Cc,Ca,[1,1,0],R);
[Bce,Bae,De]=ncdestbd(u,y,Ac,Aa,Cc,Ca,[0,1,0]);             

