clear all 
n=4;l=3;m=2;N=100;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);


% standard tests
[Be,De]=destbd(u,y,A,C);
ye=dlsim(A,Be,C,De,u);
V=vaf(y,ye);
if V==100,disp('test succeeded'),else disp('test failed');V,end
[Be,De,R]=destbd(u,y,A,C);
[Be,De,R,Phi]=destbd(u,y,A,C);

[Be,De,R,Phi]=destbd(u,y,A,C);
[Be,De,x0e,R,Phi]=destbd(u,y,A,C,[1 1 0]);
[Be,De,x0e,R,Phi]=destbd(u,y,A,C,[1 0 0]);
[Be,De,x0e,R,Phi]=destbd(u,y,A,C,[0 1 0]);
[Be,De,x0e,R,Phi]=destbd(u,y,A,C,[0 0 1]);
[Be,De,x0e,R,Phi]=destbd(u,y,A,C,[0 0 0]);
disp('test succeeded')




% test kalman filter estimate (undocumented)
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
yk = dlsim(A-K*C,[B-K*D,K],C,[D zeros(3)],[u y]);
V1=vaf(y,yk);

% estimated kalmanfilter
[Be,De,x0e,Ke,R2,Phi2]=destbd(u,y,A-K*C,C,[1 1 1 1]);
yke = dlsim(A-Ke*C,[Be-Ke*De,Ke],C,[De zeros(3)],[u y]);
V2=vaf(y,yke);
%V1 and V2 should be about equal!
if max(abs(V1-V2))<10,disp('test succeeded'),else disp('test failed');V1,V2,end

n=6;l=1;m=1;N=1000;
[A,B,C,D]=drmodel(n,l,m);
u=randn(N,m);
y=dlsim(A,B,C,D,u);
[Be,De]=destbd(u,y,A,C);
ye=dlsim(A,Be,C,De,u);
V=vaf(y,ye);
if V==100,disp('test succeeded'),else disp('test failed');V,end
