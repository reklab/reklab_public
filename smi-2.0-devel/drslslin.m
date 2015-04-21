 
function [A,B,C,D,K,hist]  =  drslslin(u,y,A,C,K,model,partype,options,slsstate) 
 % drslslin  This function performs a recursive update of a discrete time 
%           state space system, 
%                      x(k+1) = A x(k) + B u(k) 
%                      y(k)   = C x(k) + D u(k) + v(k) 
%           using the Seperable Least Squares technique. 
%           For this, only the initial estimates of A and C are needed. 
%           The same function can also be used to initialize the matrices 
%           needed to start up the recursion: Phi,Psi,dPhi, Y 
% Syntax: 
%           [A,B,C,D,hist] = drslslin(u,y,A,C) 
%           [A,B,C,D,hist] = drslslin(u,y,A,C,K,model,partype,options,slsstate) 
%           [slsstate] = drslslin(u,y,A,C,K,model,partype,'init') 
% 
% Input: 
% A,C       Initial estimate of the state space matrices A and C 
% u,y       Input 
% model     Vector [fB,fD,fx,fK] specifying whether the matrix B, D, 
%           the initial state and the kalman filter gain K 
%           should be estimated. A zero means 
%           that the parameters will not be estimated. Default is 
%           [1 0 0]. It is recommended only to estimate the initial 
%           state if the data length is short compared to the 
%           largest time constant of the system. 
% options   This vector consists of parameters that can be used to influence 
%           the optimization process. 
%           options(1) Display parameter (Default:0). 1 displays some results. 
%           options(2) is the stepsize parameter. 
%           options(3) is a forgetting factor of the nonlinear param. 
%           options(4) is a forgetting factor of the linear param. 
% 
% Outputs 
% A,B,C,D  estimated state space matrices 
% hist     history of the recursion 
%          the first column is the value of the costfunction at every step. 
%          the other colums are the estimated parameters at every step. 
% 
% See also: dslslin 
 
 
 if nargin==0 
  help drslslin 
  return 
end 
 
if nargin<4 
  error('Not enough input arguments') 
end 
 
if size(y,2)>size(y,1) 
  y = y'; 
end 
if size(u,2)>size(u,1) 
  u = u'; 
end 
 
n = length(A); 
m = size(u,2); 
N = length(u); 
 
error(abcdchk(A,[],C,[])) 
if ~(max(abs(eig(A)))<1) 
 error('Initial A must be stable.') 
end 
if rank(obsv(A,C))<length(A) 
 error('Initial (A,C) must be observable.') 
end 
n = size(A,1); 
l = size(C,1); 
 
if nargin<5 
  K = []; 
end 
 
if ~isempty(K) & size(K)~=[n,l] 
  error('Kalman matrix has wrong dimensions') 
end 
 
if nargin<6 
  model = []; 
end 
if isempty(model) 
  model = [1 1 0 0]; 
end 
if length(model)<3 
  error('Model parameter should contain 3 flags'); 
else 
  fB = model(1); 
  fD = model(2); 
  fx = model(3); 
  fK = model(4); 
end 
 
 
if nargin<7 
  partype = []; 
end 
if isempty(partype) 
  partype = 'on'; 
end 
 
if  partype=='on' 
  params.partype = partype; 
  nn = n * l; 
elseif  partype=='tr' 
  params.partype = partype; 
  nn = n * l+ 3 * n-2; 
else 
  error('You specified an unknown type of parameterization.') 
end 
 
if nargin<8 
  options = []; 
end 
if isempty(options) 
  options = [1 1 1 1]; 
end 
if length(options)==1 
  options = [options, 1 1 1]; 
elseif length(options)~=4 
  error('Options has wrong length') 
end 
 
initflag = strcmp(options,'init'); 
if initflag 
  options = [0 1 1 1]; 
end 
dispflag = options(1); 
mu = options(2); 
mu0 = mu; 
nun = options(3); 
nul = options(4); 
 
 
 
if nargin<9 
  slsstate = []; 
end 
 
nl = fx * n+fB * n * m+fD * l * m+fK * n * l; 
 
if (nargin==8)&(initflag==0) 
  init_from_slsstate_flag = 0; 
else 
  init_from_slsstate_flag = 1; 
end 
 
 
params.n = n; 
params.m = m; 
params.l = l; 
params.fx = fx; 
params.fB = fB; 
params.fD = fD; 
params.fK = fK; 
params.partype = partype; 
params.nn = nn; 
params.nl = nl; 
Inn = speye(nn); 
Il = eye(l); 
Im = eye(m); 
P1 = spalloc(nl,nl,(fx+m) * n^2); 
 
 
if initflag 
  % start initialisation mode 
    Rl_km1 = eye(nl) * 1e10; 
phi_km1 = zeros(l,nl); 
dphidt_km1 = zeros(l,nl * nn); 
u_km1 = zeros(m,1); 
y_km1 = zeros(l,1); 
Rn_km1 = speye(nn) * 1e10; 
Cl_km1 = zeros(nl,nn * nl); 
Dl_km1 = zeros(nl,nn); 
epsilon_km1 = zeros(l,1); 
thl_km1 = zeros(nl,1); 
 
thn_km1 = dss2th(A,C,params.partype); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
[B,D,x0,K]  =  destbd(u,y,A,C,model); 
theta = dss2th(A,B,C,D,x0,K,params.partype); 
thl_km1 = theta(nn+1:end); 
for k = 1:N 
  y_k = y(k,:)'; 
  u_k = u(k,:)'; 
   [A,C] = dth2ss(thn_km1,params); 
  P1 = spalloc(nl,nl,(fx+m) * n^2); 
P1(1:(fB * m+fx+fK * l) * n,1:(fB * m+fx+fK * l) * n) = kron(speye(fB * m+fx+fK * l),A); 
if ~fK 
  P2 = [zeros(l,fx * n),kron(u_km1',C),kron(u_k',Il)]; 
else 
    P2 = [zeros(l,fx * n),kron(u_km1',C),kron(y_km1',C),kron(u_k',Il)]; 
end 
phi_k = phi_km1 * P1+P2; 
 
 
 
L_k = Rl_km1 * phi_k' * inv(nul * Il+phi_k * Rl_km1 * phi_k'); 
Rl_k = (Rl_km1-L_k * phi_k * Rl_km1)/nul; 
thl_k = thl_km1 + L_k * (y_k-phi_k * thl_km1); 
 
 
 
  ye_k =  phi_k * thl_k; 
epsilon_k = y_k(:)-ye_k(:); 
 
 
  [dAdt,dCdt] = calcdadc(thn_km1,params); 
dP1dt = spalloc(nl,nn * nl,(fB * m+fx+fK * l) * n^2 * nn); 
dP2dt = spalloc(l,nn * nl,nn * l * n * m); 
w1 = [1:fx * n+fB * n * m+fK * n * l]; 
w2 = fx * n+1:fx * n+n * m; 
w3 = fx * n+n * m+1:fx * n+n * m+l * n; 
for j = 1:nn 
  dP1dt(w1,(j-1) * nl+w1) = kron(speye(fx+fB * m+fK * l),dAdt(:,(j-1) * n+1:j * n)); 
  if fB 
    dP2dt(:,(j-1) * nl+w2) = kron(u_km1',dCdt(:,(j-1) * n+1:j * n)); 
  end 
  if fK 
    dP2dt(:,(j-1) * nl+w3) = kron(y_km1',dCdt(:,(j-1) * n+1:j * n)); 
  end 
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
 
if partype=='on' &~initflag 
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
 
 
 
 
 
 
  thn_k = thn_km1; 
   Rl_km1 = Rl_k; 
phi_km1 = phi_k; 
u_km1 = u_k; 
y_km1 = y_k; 
thl_km1 = thl_k; 
dphidt_km1 = dphidt_k; 
Rn_km1 = Rn_k; 
Cl_km1 = Cl_k; 
Dl_km1 = Dl_k; 
thn_km1 = thn_k; 
Cost_k = epsilon_k' * epsilon_k; 
hist = [hist;Cost_k,thn_k',thl_k',norm(psi_k),norm(dphidt_k),  ... 
      norm(dthldt_k)]; 
if dispflag 
disp(['Step: ', num2str(k), ... 
      ' err: ', num2str(Cost_k), ... 
      ' thn: ', num2str(transpose(thn_km1(1:min(2,nn)))), ... 
      ' thl: ', num2str(transpose(thl_km1(1:min(2,nl))))]); 
end 
 
 
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
  if   init_from_slsstate_flag 
     Rl_km1 = slsstate.Rl; 
phi_km1 = slsstate.phi; 
dphidt_km1 = slsstate.dphidt; 
u_km1 = slsstate.u; 
Rn_km1 = slsstate.Rn; 
Cl_km1 = slsstate.Cl; 
Dl_km1 = slsstate.Dl; 
thl_km1 = slsstate.thl; 
thn_km1 = dss2th(A,C,params.partype); 
epsilon_km1 = zeros(l,1); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
  else 
     Rl_km1 = eye(nl) * 1e10; 
phi_km1 = zeros(l,nl); 
dphidt_km1 = zeros(l,nl * nn); 
u_km1 = zeros(m,1); 
y_km1 = zeros(l,1); 
Rn_km1 = speye(nn) * 1e10; 
Cl_km1 = zeros(nl,nn * nl); 
Dl_km1 = zeros(nl,nn); 
epsilon_km1 = zeros(l,1); 
thl_km1 = zeros(nl,1); 
 
thn_km1 = dss2th(A,C,params.partype); 
Cost_km1 = epsilon_km1' * epsilon_km1; 
hist = [Cost_km1,thn_km1',thl_km1',0,0,0]; 
 
 
 
  end 
   for k = 1:N 
  y_k = y(k,:)'; 
  u_k = u(k,:)'; 
   [A,C] = dth2ss(thn_km1,params); 
  P1 = spalloc(nl,nl,(fx+m) * n^2); 
P1(1:(fB * m+fx+fK * l) * n,1:(fB * m+fx+fK * l) * n) = kron(speye(fB * m+fx+fK * l),A); 
if ~fK 
  P2 = [zeros(l,fx * n),kron(u_km1',C),kron(u_k',Il)]; 
else 
    P2 = [zeros(l,fx * n),kron(u_km1',C),kron(y_km1',C),kron(u_k',Il)]; 
end 
phi_k = phi_km1 * P1+P2; 
 
 
 
L_k = Rl_km1 * phi_k' * inv(nul * Il+phi_k * Rl_km1 * phi_k'); 
Rl_k = (Rl_km1-L_k * phi_k * Rl_km1)/nul; 
thl_k = thl_km1 + L_k * (y_k-phi_k * thl_km1); 
 
 
 
  ye_k =  phi_k * thl_k; 
epsilon_k = y_k(:)-ye_k(:); 
 
 
  [dAdt,dCdt] = calcdadc(thn_km1,params); 
dP1dt = spalloc(nl,nn * nl,(fB * m+fx+fK * l) * n^2 * nn); 
dP2dt = spalloc(l,nn * nl,nn * l * n * m); 
w1 = [1:fx * n+fB * n * m+fK * n * l]; 
w2 = fx * n+1:fx * n+n * m; 
w3 = fx * n+n * m+1:fx * n+n * m+l * n; 
for j = 1:nn 
  dP1dt(w1,(j-1) * nl+w1) = kron(speye(fx+fB * m+fK * l),dAdt(:,(j-1) * n+1:j * n)); 
  if fB 
    dP2dt(:,(j-1) * nl+w2) = kron(u_km1',dCdt(:,(j-1) * n+1:j * n)); 
  end 
  if fK 
    dP2dt(:,(j-1) * nl+w3) = kron(y_km1',dCdt(:,(j-1) * n+1:j * n)); 
  end 
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
 
if partype=='on' &~initflag 
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
 
 
 
 
 
 
   Rl_km1 = Rl_k; 
phi_km1 = phi_k; 
u_km1 = u_k; 
y_km1 = y_k; 
thl_km1 = thl_k; 
dphidt_km1 = dphidt_k; 
Rn_km1 = Rn_k; 
Cl_km1 = Cl_k; 
Dl_km1 = Dl_k; 
thn_km1 = thn_k; 
Cost_k = epsilon_k' * epsilon_k; 
hist = [hist;Cost_k,thn_k',thl_k',norm(psi_k),norm(dphidt_k),  ... 
      norm(dthldt_k)]; 
if dispflag 
disp(['Step: ', num2str(k), ... 
      ' err: ', num2str(Cost_k), ... 
      ' thn: ', num2str(transpose(thn_km1(1:min(2,nn)))), ... 
      ' thl: ', num2str(transpose(thl_km1(1:min(2,nl))))]); 
end 
 
 
end 
 
 
 
  [A,B,C,D,x0,K] = dth2ss([thn_km1;thl_km1],params); 
end 
% add all subfunction 
 function [dAdt,dCdt] = calcdadc(thn,params) 
n = params.n; 
l = params.l; 
nn = params.nn; 
partype = params.partype; 
dAdt = zeros(n,n * nn); 
dCdt = zeros(l,n * nn); 
e = eye(l); 
if partype=='on' 
  %output normal 
  for j = 1:nn 
    Z = [zeros(l,n);eye(n)]; 
    for i = 1:n 
      si = thn(l * (n-i)+1:l * (n-i+1)); 
      ti = si' * si; 
      ri = sqrt(1-ti); 
      if ceil(j/l)==n-i+1 
 s = rem(j,l); 
 if s==0 
   s = l; 
 end 
 Ti = zeros(n+l); 
 if ti>eps 
   dSdt = (2 * thn(j) * (1-ri)/ti^2 - thn(j)/(ri * ti)) * si * si' ... 
       -(1-ri)/ti *  (e(:,s) * si'+si * e(s,:)); 
   Ti(i:l+i,i:l+i) = [-e(:,s) dSdt; -thn(j)/ri e(s,:)]; 
 else 
   Ti(i:l+i,i:l+i) = [-e(:,s) zeros(l); -thn(j) e(s,:)]; 
 end 
      else 
 Ti = eye(n+l); 
 if ti>eps 
   Ti(i:l+i,i:l+i) = [-si,e-(1-ri)/ti * si * si';ri,si']; 
 else 
   %limit ti->0 
   Ti(i:l+i,i:l+i)  =  [-si,Il-0.5 * si * si';ri,si']; 
 end 
      end 
      Z = Ti' * Z; 
    end 
    dAdt(:,(j-1) * n+1:j * n) = Z(l+1:l+n,:); 
    dCdt(:,(j-1) * n+1:j * n) = Z(1:l,:); 
  end 
elseif partype=='tr' 
  % tri-diagonal 
  dAdtj = zeros(n); 
  dCdtj = zeros(l,n); 
  for j = 1:n-1 
    dAdtj(j,j+1) = 1; 
    dAdt(:,(j-1) * n+1:j * n) = dAdtj; 
    dAdtj(j,j+1) = 0; 
  end 
  for j = 1:n 
    dAdtj(j,j) = 1; 
    dAdt(:,(n-1) * n+[(j-1) * n+1:j * n]) = dAdtj; 
    dAdtj(j,j) = 0; 
  end 
  for j = 1:n-1 
    dAdtj(j+1,j) = 1; 
    dAdt(:,(n-1) * n+n^2+[(j-1) * n+1:j * n]) = dAdtj; 
    dAdtj(j+1,j) = 0; 
  end 
  na = 3 * n-2; 
  for j = 1:n 
    for i = 1:l 
      dCdtj(i,j) = 1; 
      dCdt(:,na * n+[(j-1) * n * l+(i-1) * n+1:(j-1) * n * l+i * n]) = dCdtj; 
      dCdtj(i,j) = 0; 
    end 
  end 
end 
 
 
 function dphidt =  ncalcdphi(thn0,u,params); 
delta = 1e-6; 
nn = params.nn; 
nl = params.nl; 
l = params.l; 
dphidt = zeros(l,nl * nn); 
size(u) 
phi0 = ncalcphi(thn0,u,params); 
for j = 1:nn 
  thn1 = thn0; 
  thn1(j) = thn1(j)+delta; 
  phi1 = ncalcphi(thn1,u,params); 
  dphidt(:,(j-1) * nl+1:j * nl) = (phi1-phi0)/delta; 
end 
 
 
 function phi = ncalcphi(thn,u,params) 
[A,C] = dth2ss(thn,params); 
k = size(u,1); 
m = params.m; 
n = params.n; 
l = params.l; 
fx = params.fx; 
phi = zeros(l,n * m); 
Im = eye(m); 
Il = eye(l); 
if fx 
  CAn = C; 
  for t = 2:k 
    CAn = CAn * A; 
    phi = phi * kron(Im,A)+kron(u(t-1,:),C); 
  end 
  phi = [CAn,phi,kron(u(k,:),Il)]; 
else 
  for t = 2:k 
    phi = phi * kron(Im,A)+kron(u(t-1,:),C); 
  end 
  phi = [phi,kron(u(k,:),Il)]; 
end 
 
 
 
 
 function [dAdt,dCdt] = ncalcdadc(thn0,params) 
n = params.n; 
l = params.l; 
nn = params.nn; 
partype = params.partype; 
delta = 1e-9; 
[A0,C0] = dth2ss(thn0,params); 
dAdt = zeros(n,n * nn); 
dCdt = zeros(l,n * nn); 
for j = 1:nn 
  thn1 = thn0; 
  thn1(j) = thn1(j)+delta; 
  [A1,C1] = dth2ss(thn1,params); 
  dAdt(:,(j-1) * n+1:j * n) = (A1-A0)/delta; 
  dCdt(:,(j-1) * n+1:j * n) = (C1-C0)/delta; 
end 
 
 
 function dthldt = ncalcdthl(thn0,u,y,params) 
nn = params.nn; 
delta = 1e-6; 
Phi0 = ncalcPhi(thn0,u,params); 
Y = y';Y = Y(:); 
thl0 = Phi0 \ Y; 
for j = 1:nn 
  thn1 = thn0; 
  thn1(j) = thn1(j)+delta; 
  Phi1 = ncalcPhi(thn1,u,params); 
  thl1 = Phi1 \ Y; 
  dthldt(:,j) = (thl1-thl0)/delta; 
end 
 
 
 function psi = ncalcpsi(thn0,u,y,params); 
nn = params.nn; 
l = params.l; 
delta = 1e-6; 
psi = zeros(l,nn); 
k = size(u,1); 
Y = y';Y = Y(:); 
y_k = y(k,:)'; 
Phi0 = ncalcPhi(thn0,u,params); 
phi0 = Phi0((k-1) * l+1:k * l,:); 
thl0 = Phi0 \ Y; 
epsilon0 = y_k-phi0 * thl0; 
for j = 1:nn 
  thn1 = thn0; 
  thn1(j) = thn1(j)+delta; 
  Phi1 = ncalcPhi(thn1,u,params); 
  phi1 = Phi1((k-1) * l+1:k * l,:); 
  thl1 = Phi1 \ Y; 
  epsilon1 = y_k-phi1 * thl1; 
  psi(:,j) = (epsilon1-epsilon0)/delta; 
end 
 
 
 
 
 
 
 
 
 
  
 function Phi = ncalcPhi(thn,u,params) 
[A,C] = dth2ss(thn,params); 
k = size(u,1); 
m = params.m; 
n = params.n; 
l = params.l; 
fx = params.fx; 
phi = zeros(l,n * m); 
 
Im = eye(m); 
Il = eye(l); 
if fx 
  Phi = [C,zeros(l,n * m),kron(u(1,:),Il)]; 
  CAn = C; 
  for t = 2:k 
    CAn = CAn * A; 
    phi = phi * kron(Im,A)+kron(u(t-1,:),C); 
    Phi = [Phi;CAn,phi,kron(u(t,:),Il)]; 
  end 
else 
  Phi = [zeros(l,n * m),kron(u(1,:),Il)]; 
  for t = 2:k 
    phi = phi * kron(Im,A)+kron(u(t-1,:),C); 
    Phi = [Phi;phi,kron(u(t,:),Il)]; 
  end 
end 
 
 

