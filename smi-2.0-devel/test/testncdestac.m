%make anticausal model 
u=rand(101,1);
y=u(2:101);
u=u(1:100);
[S,R]=dordpo(u,y,4);
[Ace,Aae,Cce,Cae]=ncdestac(R,1);
eigAce=eig(Ace)
eigAae=eig(Aae)


n=4;l=3;m=2;N=101;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);
y=u(2:101);
u=u(1:100);

[S,R]=dordpo(u,y,8);
[Ace,Aae,Cce,Cae]=ncdestac(R,5);
Ace
Aae
Cce
Cae