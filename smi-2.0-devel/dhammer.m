function[A,B,C,D,NN,ze] = dhammer(u,z,n,nn);

% still to come




% checks are still to come
N = size(u,1);
m = size(u,2);
l = size(z,2);




%first build tchebychev polynomials of inputs
%shift=(max(u)+min(u))/2;
%scale=(max(u)-min(u))/2;
%u=(u-ones(N,1)*shift)./(ones(N,1)*scale);
T = [ones(size(u)),u,zeros(N,(nn-2) * m)];% zero and first order
for j = 3:nn+1
T(:,(j-1) * l+1:j * l) = 2 * u.*T(:,(j-2) * l+1:(j-1) * l)-T(:,(j-3) * l+1:(j-2) * l);
end



i = 2 * n;
[S,R] = dordpi(u,z,i);
[A,C] = destac(R,n);





%< Build chebychev polinomial from input >

% build linear output of linear system, with tchebychev polynomial as
%  input.
Zijt = zeros(N * l,n+nn * m * n);
e = eye(n);
temp = zeros(N * l,1);
% special case (zeroth order): first m columns of T are constants
% we need only one of those.
t = 1;
Tt = T(:,1);
for i = 1:n
xij = ltitr(A,e(:,i),Tt,zeros(n,1));
yij = xij * C';
temp(:) = yij;
Zijt(:,i) = temp;
end
for t = 2:nn+1
Tt = T(:,(t-1) * m+1:t * m);
for j = 1:m
for i = 1:n
xijt = ltitr(A,e(:,i),Tt(:,j),zeros(n,1));
yijt = xijt * C';
temp(:) = yijt;
Zijt(:,n+(t-2) * m * n+(j-1) * n+i) = temp;
end
end
end

%Uijt
nt = 1+nn * m;
Uijt = zeros(N * l:nt * l);
for j = 1:l,
Uijt((j-1) * N+1:j * N,(j-1) * nt+1:j * nt) = T(:,m:m+nn * m);
end

Zijt(1:5,:)
condzijt = cond(Zijt)
sizezijt = size(Zijt)
phi = [Zijt,Uijt];





NN = pinv(phi) * z(:);
BNN = zeros(n,1+nn * m);
BNN(:) = NN(1:n+nn * n * m);
DNNT = zeros(1+nn * m,l);
DNNT(:) = NN(n+nn * n * m+1:n+nn * n * m+l+nn * l * m);DNN = DNNT';


B = BNN;
D = DNN;
ze = phi * NN;
ze = reshape(ze,N,l);






