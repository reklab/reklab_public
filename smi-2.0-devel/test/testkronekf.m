A=randn(5);
B=randn(5);
[AA,BB,Q,Z,na,nb]=kronekf(A,B);
AA-Q*A*Z
BB-Q*B*Z
eig(inv(BB)*AA)
eig(inv(B)*A)


A=randn(5);
B=eye(5);
[AA,BB,Q,Z,na,nb]=kronekf(A,B);
eigA=eig(A)
inveigBB=1./eig(BB)
eigAA=eig(AA)

[AA,BB,Q,Z,na,nb]=kronekf(A,B,2,'circle');
[AA,BB,Q,Z,na,nb]=kronekf(A,B,1,'circle');

[AA,BB,Q,Z,na,nb]=kronekf(A,B,0,'halfplane');

eig(A,B)
eig(AA,BB)
eig(AA)
eig(BB)





