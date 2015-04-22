n=4;l=3;m=2;N=100;i=6
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);

beta=0.9
[S,R]=dordpo(u,y,i);
Ai=A+0.2*randn(size(A));
Bi=B+0.2*randn(size(B));
Ci=C+0.2*randn(size(C));
Di=D+0.2*randn(size(D));
yi=dlsim(Ai,Bi,Ci,Di,u);


[Ao,Bo,Co,Do]= drpo(u,y,i,n,beta,Ai,Bi,Ci,Di,R); 
 eig(Ao(:,137:140))  
eig(A)

Ae=Ao(:,137:140);
Be=Bo(:,69:70);
Ce=Co(:,137:140);
De=Do(:,69:70);
ye=dlsim(Ae,Be,Ce,De,u);
vaf(y,yi)
vaf(y,ye)
