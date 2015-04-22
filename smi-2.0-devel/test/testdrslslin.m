n=4;l=3;m=2;N=200;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


Ai=A+0.2*randn(size(A));
Ci=C+0.2*randn(size(C));
options=1;
partype='on';
[Ae,Be,Ce,De,hist] = drslslin(u,y,Ai,Ci,[],[],[],0);
ye1=dlsim(Ai,B,Ci,D,u);
V1=vaf(y,ye1)
ye2=dlsim(Ae,Be,Ce,De,u);
V2=vaf(y,ye2)


if 0
% test kalman filter estimate
% this is probably one step to far for this function;-(
% the function runs, but the answers are crappy.
% I haven't found a way to test it on a trivial example,
% due to the noise, it is always an approximation...
N=1000;N1=100;a=0.1;
u=rand(N,m);
Q=a*eye(n);
R=a*eye(l);
G=a*eye(n);
v=rand(N,l);
w=rand(N,n);
y = dlsim(A,[B G],C,[D zeros(l,n)],[u, w]);

%theoretical kalman filter
L = dlqe(A,G,C,Q,R);
K = A*L;
yk = dlsim(A-K*C,[B-K*D,K],C,[D zeros(3)],[u y]);
V1=vaf(y,yk)
alpha=0;
Ai=A+alpha*randn(size(A));
Ci=C+alpha*randn(size(C));
Ki=K+alpha*randn(size(K));
model=[1 1 0 1];
options=0;
partype='on';

[slsstate] = drslslin(u(1:N1,:),y(1:N1,:),A,C,K,model,partype,'init')
[Ae,Be,Ce,De,Ke,hist] = drslslin(u(N1+1:N,:),y(N1+1:N,:),Ai,Ci,Ki,model,[],0);

ye2=dlsim(Ae,Be,Ce,De,u);
V2=vaf(y,ye2)

end