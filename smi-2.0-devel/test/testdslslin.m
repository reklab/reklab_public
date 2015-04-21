clear all
n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


[Ae,Be,Ce,De]=dslslin(u,y,A,C);
ye=dlsim(Ae,Be,Ce,De,u);
V=vaf(y,ye);
if V==100,disp('test succeeded'),else disp('test failed');V,end


[Ae,Be,Ce,De,x0e]=dslslin(u,y,A,C,[],[1 1 0 0],'on');
[Ae,Be,Ce,De,x0e]=dslslin(u,y,A,C,[],[1 0 0 0],'on');
[Ae,Be,Ce,De,x0e]=dslslin(u,y,A,C,[],[0 1 0 0],'on');
[Ae,Be,Ce,De,x0e]=dslslin(u,y,A,C,[],[0 0 1 0],'on');
[Ae,Be,Ce,De,x0e]=dslslin(u,y,A,C,[],[0 0 0 0],'on');
disp('test succeeded')

% test kalman filter estimate
N=1000;
u=rand(N,m);
Q=eye(n);
R=eye(l);
G=eye(n);
v=rand(N,l);
w=rand(N,n);
y = dlsim(A,[B G],C,[D zeros(l,n)],[u, w]);

%theoretical kalman filter
L = dlqe(A,G,C,Q,R);
K = A*L;
yk = dlsim(A-K*C,[B-K*D,K],C,[D zeros(l)],[u y]);
V1=vaf(y,yk)

% Start value of dslslin
%[Ba,Da,xa,Ka]=destbd(u,y,A-K*C,C,[1 1 0 1]);
%yka = dlsim(A-K*C,[Ba-Ka*Da,Ka],C,[Da zeros(l)],[u y]);
%V1a=vaf(y,yka)


[Ae,Be,Ce,De,x0e,Ke] = dslslin(u,y,A,C,K,[1 1 0 1],[],1);
yke = dlsim(Ae-Ke*Ce,[Be-Ke*De,Ke],Ce,[De zeros(l)],[u y]);
V2=vaf(y,yke)








