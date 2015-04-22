n=4;l=3;m=2;N=1000;i=6;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


% test kalman filter estimate
N=1000;
u=rand(N,m);
Q=eye(n);
R=eye(l);
G=eye(n);
v=rand(N,l);
w=rand(N,n);
y = dlsim(A,[B G],C,[D zeros(l,n)],[u, 0.5*w])+0.5*v;

%theoretical kalman filter
L = dlqe(A,G,C,Q,R);
K = A*L;
yk = dlsim(A-K*C,[B-K*D,K],C,[D zeros(3)],[u y]);
V1=vaf(y,yk)

% estimated kalmanfilter
[Sn,R]=dordpo(u,y,i); 
[Ae,Ce]=destac(R,n);
[Be,De]=destbd(u,y,Ae,Ce);
Ke=destk(Ae,Ce,R);

yke = dlsim(Ae-Ke*Ce,[Be-Ke*De,Ke],Ce,[De zeros(3)],[u y]);
V1=vaf(y,yke)
%V1 and V2 should be about equal!




