function bl  = nlident (bl, z, varargin)
% Identify a lnbl

if nargin > 2
   set (bl,varargin{:});
end

if isa(z,'nldat') | isa(z,'double')
  if isa(z,'nldat')
    Ts=z.domainIncr;
  else  
    subsys = get(bl,'elements');
    f1 = subsys{1};
    Ts = f1.domainIncr;
    z = nldat(z,'domainIncr',Ts);
  end   

assign(bl.parameterSet);
subsys = bl.elements;
  i = subsys{1,1};  % IRF
  p = subsys{1,2};  % Polynomial
  x=z(:,1);
  y=z(:,2);
  hlen = nLags;
  set(i,'nLags',hlen,'orderSelectMode',orderSelectMode);
  switch lower(initMethod)
    case 'correl'
      phi = cor('kernOrder',1,'nLags',hlen);
      phi = nlident(phi,z);
      h_est = i;
      set(h_est,'dataSet',double(phi));    
    case 'fil'
      h_est = nlident(i, z);
    case 'slice2'
      h_est = cor_slice(i,z,2);
    case 'slice3' 
      h_est = cor_slice(i,z,3);    
    case 'eigen'  
      phi = cor;
      set(phi,'kernOrder',2,'nLags',hlen);
      phi = nlident(phi,z);
      [U,S,V] = svd(double(phi));
      h_est = i;
      set(h_est,'dataSet',U(:,1),'domainIncr',z.domainIncr);    
    case 'gen_eigen'
      h_est = wiener_2(z,i);
    otherwise
      error(['unrecognized initialization method': initMethod]);
  end


  x_est = nlsim (h_est,x);
  z1 = cat(2,x_est,y);
  set (p,'polyOrderMax',polyOrderMax,'polyOrderMax',Q,'polyOrderSelectMode',mode);
  m_est = nlident(p,z1);
  yp = nlsim(m_est,x_est);
  vf = vaf(y,yp);
  set (h_est,'Comment','Linear element');
  set (m_est,'Comment','Static NonLinear Element');
  set (bl,'Elements', { h_est m_est});
  


  switch Method
    case 'bussgang'
      set(bl,'Comment','LN model identified using Busgang''s theorm');
    case 'hk'
      bl = hk_ident(bl,z);
    case 'phk'
      bl = phk_ident(bl,z);
    case 'lm'
      bl = lm_ident(bl,z);
    otherwise
      error('unrecognized identification method');
  end

else
  error('conversions to models of class lnbl not yet implemented');
end


return


function h = cor_slice(h,z,order)

assign(h.Parameters); 
ud = double(z(:,1));
yd = double(z(:,2));
N = length(ud);
uny = z;

switch order
  case 2    
    lag = floor(hlen*rand(1));
    udel = [zeros(lag,1);ud(1:N-lag)];
    set(uny,'Data',[ud udel.*yd]);
  case 3
    lag1 = floor(hlen*rand(1));
    ud1 = [zeros(lag1,1);ud(1:N-lag1)];
    lag2 = floor(hlen*rand(1));
    ud2 = [zeros(lag2,1);ud(1:N-lag2)];
    set(uny,'Data',[ud ud1.*ud1.*yd]);
  otherwise
    error('unsupported slice order');
end

phi = cor;
set(phi,'Order',1,'nLags',hlen);
phi = nlident(phi,uny);
hd = phi.Data;


gain = std(yd);
switch order
  case 2    
    hd(lag+1) = hd(lag+1) + randn(1)*gain;
  case 3
    hd(lag1+1) = hd(lag1+1) + randn(1)*gain;
    hd(lag2+1) = hd(lag2+1) + randn(1)*gain;
  otherwise
    error('unsupported slice order');
end

set(h,'data',hd);

return