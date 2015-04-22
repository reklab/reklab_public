 
function [A,B,C,D,hist]  =  crslslin(A,C,K,u,y,model,options,slsstate) 
 % crslslin  This function performs a recursive update of a discrete time 
%           state space system, 
%                      x(k+1) = A x(k) + B u(k) 
%                      y(k)   = C x(k) + D u(k) + v(k 
%           using the Seperable Least Squares technique. 
%           For this, only the initial estimates of A and C are needed. 
%           The optimization is a 
%           This function can also be used to initialize the matrices 
%           needed to start up the recursion: Phi,Psi,dPhi, Y 
% Syntax: 
%           [A,B,C,D,hist] = drslslin(A,C,u,y,model,options) 
%           [A,B,C,D,hist] = drslslin(A,C,u,y,model,options,slsstate) 
%           [slsstate] = drslslin(A,C,u,y,model,'init') 
% 
% Input: 
% A,C       Initial estimate of the state space matrices A and C 
% u,y       Input 
% model     Vector [fD,fx,fK] specifying whether the matrix D, 
%           the initial state and the kalman filter gain K 
%           should be estimated. A zero means 
%           that the parameters will not be estimated. Default is 
%           [1 0 0]. It is recommended only to estimate the initial 
%           state if the data length is short compared to the 
%           largest time constant of the system. 
%  options   This vector consists of parameters that can be used to influence 
%            the optimization process. 
%            options(1) is the stepsize parameter. 
%            options(2) is a forgetting factor of the nonlinear param. 
%            options(3) is a forgetting factor of the linear param. 
% 
% Outputs 
% A,B,C,D  estimated state space matrices 
% hist     history of the recursion 
%          the first column is the value of the costfunction at every step. 
%          the other colums are the estimated parameters at every step. 
 
 
% function [A,B,C,D,hist] = drslslin(A,C,u,y,model,options) 
% function [A,B,C,D,hist] = drslslin(A,C,u,y,model,options,slsstate) 
% function [Psi,Phi,dPhiTdt,Y] = drslslin(A,C,u,y,model,'init') 
 
 
 
 if nargin==0 
  help crslslin 
  return 
end 
if nargin<3 
  error('Not enough input arguments') 
end 
 
error(abcdchk(A,[],C,[])) 
n = size(A,1); 
l = size(C,1); 
 
% decide if Kalman filter must be estimated 
if all(size(K)==[n,l]) 
  fK = 1; 
else 
  fK = 0; 
end 
 
if fK 
  if nargin<8 
    slsstate = []; 
  end 
  if nargin<7 
    options = []; 
  end 
  if nargin<6 
    model = []; 
  end 
else 
  if nargin<7 
    slsstate = []; 
  else 
    slsstate = options; 
  end 
 
  if nargin<6 
    options = []; 
  else 
    options = model; 
  end 
  if nargin<5 
    model = []; 
  else 
    model = y; 
  end 
  y = u; 
  u = K; 
end 
 
if size(y,2)>size(y,1) 
  y = y'; 
end 
if size(u,2)>size(u,1) 
  u = u'; 
end 
 
n = length(A); 
 
l = size(y,2); 
m = size(u,2); 
N = length(u); 
 
if nargin<5 
  model = []; 
end 
if isempty(model) 
  model = [1,0,0,0]; 
elseif length(model)<4 
  error('Model paramgeter should contain 4 flags'); 
else 
  fD = ~(model(1)==0); 
  fx = ~(model(2)==0); 
  fK = ~(model(3)==0); 
  param = model(4); 
  model = [fD,fx,fK,param]; 
end 
 
if nargin<6 
  options = []; 
end 
 
% test if options is the string 'init' 
if length(options)==4 
  init_flag = (options=='init'); 
  options = []; 
else 
  init_flag = 0; 
end 
 
% test normal use of options 
if length(options)==3 
  mu = options(1); 
  mu0 = mu; 
  if isempty(mu); 
    mu = 1; 
  end 
  nun = options(2); 
  if isempty(nun); 
    nun = 0; 
  end 
  nul = options(2); 
  if isempty(nul); 
    nul = 0; 
  end 
elseif isempty(options) 
  mu = 1; 
  nun = 1; 
  nul = 1; 
else 
  error('options has wrong size'); 
end 
 
if param==1 
  disp('tridiagonal model') 
else 
  disp('outputnormal model') 
end 
 
if param==1 
  nn = n * l+ 3 * n-2; 
else 
  nn = n * l; 
end 
nl = fx * n+n * m+fD * l * m+fK * n * l; 
 
 
if (nargin<7)&(init_flag==0) 
  init_from_slsstate_flag = 0; 
else 
  init_from_slsstate_flag = 1; 
end 
 
params.n = n; 
params.m = m; 
params.l = l; 
params.fx = fx; 
params.fD = fD; 
params.fK = fK; 
params.nn = nn; 
params.nl = nl; 
 
Inn = speye(nn); 
Il = eye(l); 
Im = eye(m); 
P1 = spalloc(nl,nl,(fx+m) * n^2); 
 
 
if init_flag 
  % start initialisation mode 
    Rl_km1 = eye(nl) * 1e10; 
phi_km1 = zeros(l,nl); 
dphidt_km1 = zeros(l,nl * nn); 
u_km1 = zeros(m,1); 
Rn_km1 = speye(nn) * 1e10; 
Cl_km1 = zeros(nl,nn * nl); 
Dl_km1 = zeros(nl,nn); 
epsilon_km1 = zeros(l,1); 
thl_km1 = zeros(nl,1); 
thn_km1 = dac2thn(A,C,params); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
for k = 1:N 
  y_k = y(k,:)'; 
  u_k = u(k,:)'; 
   [A,C] = dthn2ac(thn_km1,params); 
  P1 = spalloc(nl,nl,(fx+m) * n^2); 
P1(1:(m+fx) * n,1:(m+fx) * n) = kron(speye(m+fx),A); 
P2 = [zeros(l,fx * n),kron(u_km1',C),kron(u_k',Il)]; 
phi_k = phi_km1 * P1+P2; 
 
 
 
L_k = Rl_km1 * phi_k' * inv(nul * Il+phi_k * Rl_km1 * phi_k'); 
Rl_k = (Rl_km1-L_k * phi_k * Rl_km1)/nul; 
thl_k = thl_km1 + L_k * (y_k-phi_k * thl_km1); 
 
 
 
  ye_k =  phi_k * thl_k; 
epsilon_k = y_k(:)-ye_k(:); 
 
 
 
  [dAdt,dCdt] = calcdadc(thn_km1,params); 
dP1dt = spalloc(nl,nn * nl,(m+fx) * n^2 * nn); 
dP2dt = spalloc(l,nn * nl,nn * l * n * m); 
w1 = [1:fx * n+n * m]; 
w2 = fx * n+1:fx * n+n * m; 
for j = 1:nn 
  dP1dt(w1,(j-1) * nl+w1) = kron(speye(fx+m),dAdt(:,(j-1) * n+1:j * n)); 
  dP2dt(:,(j-1) * nl+w2) = kron(u_km1',dCdt(:,(j-1) * n+1:j * n)); 
end 
dphidt_k = dphidt_km1 * kron(Inn,P1)+phi_km1 * dP1dt+dP2dt; 
 
 
 
 dphiTdt_k = zeros(nl,l * nn); 
for j = 1:nn 
  dphiTdt_k(:,(j-1) * l+1:j * l) = dphidt_k(:,(j-1) * nl+1:j * nl)'; 
end 
Cl_k = Cl_km1+phi_k' * dphidt_k+dphiTdt_k * kron(Inn,phi_k); 
Dl_k = Dl_km1+dphiTdt_k * kron(Inn,y_k); 
dthldt_k = Rl_k * (Dl_k-Cl_k * kron(Inn,thl_k)); 
 
 
 
psi_k = -dphidt_k * kron(Inn,thl_k)-phi_k * dthldt_k; 
 
 
 
Rn_km1 = Rn_km1/nun; 
Rn_k = Rn_km1-Rn_km1 * psi_k' * inv(Il * nun+psi_k * Rn_km1 * psi_k') * psi_k * Rn_km1; 
thn_k = thn_km1-mu * Rn_k * psi_k' * epsilon_k; 
if param==0 &~init_flag 
  % check new thetan 
  i = 1; 
  flag_badthn = 0; 
  while (i<=n) & (flag_badthn==0) 
    si = thn_k(l * (n-i)+1:l * (n-i+1)); 
    ti = si' * si; 
    if (ti>1) 
      disp(['Step ', num2str(k),': warning: new thn invalid']) 
      thn_k = thn_km1; 
      flag_badthn = 1; 
    end 
    i = i+1; 
  end 
  if flag_badthn 
    mu = mu * .99; 
  else 
    mu = mu * 1.01; 
  end 
  mu = min(mu0,mu); 
end 
 
 
 
 
 
if 0 
  % test 
  [ndAdt,ndCdt] = ncalcdadc(thn_km1,params) 
  nphi  =  ncalcphi(thn_km1,u(1:k,:),params); 
  %ndphidt=ncalcdphi(thn_km1,u(1:k,:),params); 
  %[ndAdt,ndCdt]=ncalcdadc(thn_km1,params); 
  %ndthldt=ncalcdthl(thn_km1,u(1:k,:),y(1:k,:),params); 
  %npsi=ncalcpsi(thn_km1,u(1:k,:),y(1:k,:),params); 
  dCdt 
  ndCdt 
  pause 
end 
 
 
  thn_k = thn_km1; 
 
   Rl_km1 = Rl_k; 
phi_km1 = phi_k; 
u_km1 = u_k; 
thl_km1 = thl_k; 
dphidt_km1 = dphidt_k; 
Rn_km1 = Rn_k; 
Cl_km1 = Cl_k; 
Dl_km1 = Dl_k; 
thn_km1 = thn_k; 
Cost_k = epsilon_k' * epsilon_k; 
hist = [hist;Cost_k,thn_k',thl_k',norm(psi_k),norm(dphidt_k),norm(dthldt_k)]; 
disp(['Step: ', num2str(k), ... 
      ' err: ', num2str(Cost_k), ... 
      ' thn: ', num2str(transpose(thn_km1(1:min(2,nn)))), ... 
      ' thl: ', num2str(transpose(thl_km1(1:min(2,nl))))]); 
 
 
 
end 
slsstate = struct('Rl',Rl_km1, ... 
    'phi',phi_km1, ... 
    'dphidt',dphidt_km1, ... 
    'u',u_km1, ... 
    'Rn',Rn_km1, ... 
    'Cl',Cl_km1, ... 
    'Dl',Dl_km1, ... 
    'thl',thl_km1); 
A = slsstate; 
B = hist; 
 
 
else 
  % start optimization mode 
  if init_from_slsstate_flag 
     Rl_km1 = slsstate.Rl; 
phi_km1 = slsstate.phi; 
dphidt_km1 = slsstate.dphidt; 
u_km1 = slsstate.u; 
Rn_km1 = slsstate.Rn; 
Cl_km1 = slsstate.Cl; 
Dl_km1 = slsstate.Dl; 
thl_km1 = slsstate.thl; 
thn_km1 = dac2thn(A,C,params); 
epsilon_km1 = zeros(l,1); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
  else 
     Rl_km1 = eye(nl) * 1e10; 
phi_km1 = zeros(l,nl); 
dphidt_km1 = zeros(l,nl * nn); 
u_km1 = zeros(m,1); 
Rn_km1 = speye(nn) * 1e10; 
Cl_km1 = zeros(nl,nn * nl); 
Dl_km1 = zeros(nl,nn); 
epsilon_km1 = zeros(l,1); 
thl_km1 = zeros(nl,1); 
thn_km1 = dac2thn(A,C,params); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
  end 
   for k = 1:N 
  y_k = y(k,:)'; 
  u_k = u(k,:)'; 
   [A,C] = dthn2ac(thn_km1,params); 
  P1 = spalloc(nl,nl,(fx+m) * n^2); 
P1(1:(m+fx) * n,1:(m+fx) * n) = kron(speye(m+fx),A); 
P2 = [zeros(l,fx * n),kron(u_km1',C),kron(u_k',Il)]; 
phi_k = phi_km1 * P1+P2; 
 
 
 
L_k = Rl_km1 * phi_k' * inv(nul * Il+phi_k * Rl_km1 * phi_k'); 
Rl_k = (Rl_km1-L_k * phi_k * Rl_km1)/nul; 
thl_k = thl_km1 + L_k * (y_k-phi_k * thl_km1); 
 
 
 
  ye_k =  phi_k * thl_k; 
epsilon_k = y_k(:)-ye_k(:); 
 
 
 
  [dAdt,dCdt] = calcdadc(thn_km1,params); 
dP1dt = spalloc(nl,nn * nl,(m+fx) * n^2 * nn); 
dP2dt = spalloc(l,nn * nl,nn * l * n * m); 
w1 = [1:fx * n+n * m]; 
w2 = fx * n+1:fx * n+n * m; 
for j = 1:nn 
  dP1dt(w1,(j-1) * nl+w1) = kron(speye(fx+m),dAdt(:,(j-1) * n+1:j * n)); 
  dP2dt(:,(j-1) * nl+w2) = kron(u_km1',dCdt(:,(j-1) * n+1:j * n)); 
end 
dphidt_k = dphidt_km1 * kron(Inn,P1)+phi_km1 * dP1dt+dP2dt; 
 
 
 
 dphiTdt_k = zeros(nl,l * nn); 
for j = 1:nn 
  dphiTdt_k(:,(j-1) * l+1:j * l) = dphidt_k(:,(j-1) * nl+1:j * nl)'; 
end 
Cl_k = Cl_km1+phi_k' * dphidt_k+dphiTdt_k * kron(Inn,phi_k); 
Dl_k = Dl_km1+dphiTdt_k * kron(Inn,y_k); 
dthldt_k = Rl_k * (Dl_k-Cl_k * kron(Inn,thl_k)); 
 
 
 
psi_k = -dphidt_k * kron(Inn,thl_k)-phi_k * dthldt_k; 
 
 
 
Rn_km1 = Rn_km1/nun; 
Rn_k = Rn_km1-Rn_km1 * psi_k' * inv(Il * nun+psi_k * Rn_km1 * psi_k') * psi_k * Rn_km1; 
thn_k = thn_km1-mu * Rn_k * psi_k' * epsilon_k; 
if param==0 &~init_flag 
  % check new thetan 
  i = 1; 
  flag_badthn = 0; 
  while (i<=n) & (flag_badthn==0) 
    si = thn_k(l * (n-i)+1:l * (n-i+1)); 
    ti = si' * si; 
    if (ti>1) 
      disp(['Step ', num2str(k),': warning: new thn invalid']) 
      thn_k = thn_km1; 
      flag_badthn = 1; 
    end 
    i = i+1; 
  end 
  if flag_badthn 
    mu = mu * .99; 
  else 
    mu = mu * 1.01; 
  end 
  mu = min(mu0,mu); 
end 
 
 
 
 
 
if 0 
  % test 
  [ndAdt,ndCdt] = ncalcdadc(thn_km1,params) 
  nphi  =  ncalcphi(thn_km1,u(1:k,:),params); 
  %ndphidt=ncalcdphi(thn_km1,u(1:k,:),params); 
  %[ndAdt,ndCdt]=ncalcdadc(thn_km1,params); 
  %ndthldt=ncalcdthl(thn_km1,u(1:k,:),y(1:k,:),params); 
  %npsi=ncalcpsi(thn_km1,u(1:k,:),y(1:k,:),params); 
  dCdt 
  ndCdt 
  pause 
end 
 
 
   Rl_km1 = Rl_k; 
phi_km1 = phi_k; 
u_km1 = u_k; 
thl_km1 = thl_k; 
dphidt_km1 = dphidt_k; 
Rn_km1 = Rn_k; 
Cl_km1 = Cl_k; 
Dl_km1 = Dl_k; 
thn_km1 = thn_k; 
Cost_k = epsilon_k' * epsilon_k; 
hist = [hist;Cost_k,thn_k',thl_k',norm(psi_k),norm(dphidt_k),norm(dthldt_k)]; 
disp(['Step: ', num2str(k), ... 
      ' err: ', num2str(Cost_k), ... 
      ' thn: ', num2str(transpose(thn_km1(1:min(2,nn)))), ... 
      ' thl: ', num2str(transpose(thl_km1(1:min(2,nl))))]); 
 
 
 
end 
 
 
 
  [A,B,C,D,x0,K] = cth2ss([thn_km1;thl_km1],params) 
end 
% add all subfunction 

