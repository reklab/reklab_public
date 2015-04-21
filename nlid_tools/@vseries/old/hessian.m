function H = hessian(vs,z,varargin);
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
y = dz(:,2);

switch lower(method)
  case 'foa'
    
    hlen = get(vs,'NLags');
    Qmax = get(vs,'OrderMax');
    [H,G] = foa_init(u,y,hlen,Qmax);
    
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
    H = hessian(static,[u_basis y]);
    
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







