function J = jacobian(vs,z,varargin);
% computes hessian for least squares estimate of a volterra series


% Copyright 2004, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


% make sure that the second argumment is either nldat or double,
% and contains enough colums for input(s) and output


if ~(isa (z,'double') | isa(z,'nldat'))
  error (' Second argument must be of class double, or nldat');
end


method = get(vs,'method');
dz = double(z);
u = dz(:,1);


switch lower(method)
  case 'foa'
    
    hlen = get(vs,'NLags');
    Qmax = get(vs,'OrderMax');
    J = foa_jacobian(u,hlen,Qmax);   
 
    
  case 'laguerre'
    
    
    Qmax = get(vs,'OrderMax');
    hlen = get(vs,'NLags');
    alpha = get(vs,'alpha');
    num_filt = get(vs,'NumFilt');
    delay = get(vs,'Delay');

    u_basis = laguerre_basis(u,num_filt,alpha,delay);

    static = polynom;
    set(static,'Order',Qmax,'OrderMax',Qmax,...
	'NInputs',num_filt,'type','power');
    J = jacobian(static,[u_basis]);
    
  otherwise
    error('unsupported identification method');
end



function basis = laguerre_basis(x,num_filt,alpha,delay)

data_len = length(x);
x = x(:);
x = [zeros(delay,1);x(1:data_len - delay)];
basis = zeros(data_len,num_filt);
alpha_srt = sqrt(alpha);
laga = [alpha_srt,-1];
lagb = [1, -alpha_srt];

basis(:,1) = filter(sqrt(1-alpha),[1, -alpha_srt],x);
for i = 2:num_filt
  basis(:,i) = filter(laga,lagb,basis(:,i-1));
end


function J = foa_jacobian(x,nlags,Q);

N = length(x);
npar = prod(nlags:nlags+Q)/(nlags*factorial(Q));
J = zeros(N,npar);

% zero order kernel
J(:,1) = ones(N,1);

% first order kernel
if Q > 0
  J(:,2) = x;
  for i = 2:nlags
    J(:,i+1) = [0;J(1:N-1,i)];
  end
end

% second order kernel
if Q > 1
  offset = nlags+1;
  J(:,offset+1:offset+nlags) = J(:,2:nlags+1).*(x*ones(1,nlags));
  oldfirst = nlags+2;
  oldlast = 2*nlags+1;

  for i = 2:nlags
    first = oldlast+1;
    last = first+nlags-i;
    J(:,first:last) = [zeros(1,nlags-i+1);J(1:N-1,oldfirst:oldlast-1)];
    oldfirst = first;
    oldlast = last;
  end
end

    






