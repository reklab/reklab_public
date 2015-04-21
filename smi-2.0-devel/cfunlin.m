 function epsilon = cfunlin(thn,u,y,Ts,params) 
N = size(u,1); 
fK = params.fK; 
fB = params.fB; 
fD = params.fD; 
fx = params.fx; 
[A,C] = cth2ss(thn,params); 
[B,D,x0,K] = cestbd(u,y,Ts,A,C,[fB,fD,fx,fK]); 
if fK 
  ye = tlsim(A,B-fK * K * D,C,D,u,(0:N-1)' * Ts); 
else 
  ye = tlsim(A,B,C,D,u,(0:N-1)' * Ts); 
end 
epsilon = y(:)-ye(:); 
epsilon = epsilon/N; 
 
 
 
 

