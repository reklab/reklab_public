 function epsilon = dfunlin(thn,u,y,params) 
% dfunlin  This function implements the cost-fuction for dslslin 
%          It is not meant for stand-alone use. 
% 
% Syntax: 
%          epsilon=dfunlin(thn,u,y,params) 
% 
% Inputs: 
%  thn     Parameter vector describing the system matrices A and C. 
%  u,y     The input and output data of the system to be optimized. 
%  params   A structure that contains the dimension parameters of 
%          the system, such as the order, the number of inputs, 
%          whether D, x0 or K is present in the model, etc. 
% 
% Outputs: 
%  epsilon Output of the cost-function, which is the square of the error 
%          between the output and the estimated output,  divided by the 
%          number of samples. 
% 
% Seealso: dslslin dslswie dfunwie 
if nargin==0 
  help dfunlin 
  return 
end 
N = size(u,1); 
l = params.l; 
n = params.n; 
m = params.m; 
fD = params.fD; 
fx = params.fx; 
fK = params.fK; 
 
%make sure theta is valid 
thn0 = thn; 
invalid_theta_flag = 0; 
for i = 1:n 
  si = thn(l * (n-i)+1:l * (n-i+1)); 
  ti = si' * si; 
  if ti>1 
    % warning('invalid theta found in dfunlin') 
    thn(l * (n-i)+1:l * (n-i+1)) = 0.99 * si/sqrt(ti); 
    invalid_theta_flag = 1; 
  end 
end 
[A,C] = dth2ss(thn,params); 
fB = params.fB; 
fD = params.fD; 
fx = params.fx; 
fK = params.fK; 
[B,D,x0,K] = destbd(u,y,A,C,[fB,fD,1,fK]); 
if ~fB 
  B = zeros(n,m); 
end 
if ~fD 
  D = zeros(l,m); 
end 
if fK 
  ye  =  dlsim(A,[B-K * D,K],C,[D zeros(l)],[u y]); 
else 
  ye = dlsim(A,B,C,D,u); 
end 
 
epsilon = y(:)-ye(:); 
epsilon = epsilon/N; 
if invalid_theta_flag 
  epsilon = epsilon +mean(std(y)) * norm(thn-thn0); 
end 
 
 

