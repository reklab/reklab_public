 function epsilon = dfunwie(thn,u,z,params) 
 
% dfunwie  This function implements the cost-fuction for dslswie 
%          It is not meant for stand-alone use. 
% 
% Syntax: 
%          epsilon=dfun wie(thn,u,z,params) 
% 
% Inputs: 
%  thn     Parameter vector describing the system matrices A, B, C and D. 
%  u,z    The input and output data of the system to be optimized. 
%  params  A structure that contains the dimension parameters of 
%          the system, such as the order, the number of inputs, 
%          whether D or x0 is present in the model, etc. 
% 
% Outputs: 
%  epsilon Output of the cost-function, which is the square of the error 
%          between the output and the estimated output,  divided by the 
%          number of samples. 
% 
% Seealso: dslslin dslswie dfunlin 
if nargin==0 
  help dfunwie 
  return 
end 
 
N = size(u,1); 
nn = params.nn; 
fB = params.fB; 
fD = params.fD; 
fx = params.fx; 
n = params.n; 
m = params.m; 
l = params.l; 
 
%make theta valid 
thn0 = thn; 
invalid_theta_flag = 0; 
for i = 1:n 
  si = thn(l * (n-i)+1:l * (n-i+1)); 
  ti = si' * si; 
  if ti>1 
    warning('invalid theta found in dfunlin') 
    thn(l * (n-i)+1:l * (n-i+1)) = 0.99 * si/sqrt(ti); 
    invalid_theta_flag = 1; 
  end 
end 
 
[A,B,C,D,x0] = dth2ss(thn,params); 
if ~fB 
  B = zeros(n,m); 
end 
if ~fD 
  D = zeros(l,m); 
end 
ye = dlsim(A,B,C,D,u); 
 
[thl,ze] = chebest(ye,z,nn); 
epsilon = z(:)-ze(:); 
epsilon = epsilon/N; 
 
if invalid_theta_flag 
  epsilon = epsilon +mean(std(z)) * norm(thn-thn0); 
end 
 
 
 

