function pc = pcm_nlls(pc,z);
% nonlinear least squares optimization methods to fit parallel cascade
% models.  Normally called from nlident.m


% 1 extract data and tuning parameters.

u = z(:,1);
y = z(:,2);
N = length(y);
Ts = get(u,'domainincr');

P = get(pc,'parameters');
assign(P);
% creates variables containing parameters:
% Method, NPaths, NLags, OrderMax, MaxIts, Threshold, 
% accel, decel, delta, initialization  


if isnan(NPaths)
% default setting is NPaths = NaN.  This works for parallel cascade
% algorithms, as NPaths can be determined on the fly. 
% Since NLLS algorithms require a fixed structure, we will have to set
% a default.  This computation gives at least 10 data points per parameter.
  NPaths = floor(N/(10*(NLags+OrderMax)));
end

% initialize with random IRF weights and polynomial coefficients, if
% requested. 
if strcmp(initialization,'random')
  pc = init_random(pc,Ts);
end


y_est = nlsim(pc,u);
err = y - y_est;
vaf_old = double(vaf(y,y_est));

count = 1;
update_jacobian = 1;
pc_new = pc;
 
disp(['Iteration       VAF   Improvement']);

while count <= MaxIts
  if update_jacobian
    J = Jacobian(pc,u);
  end
  switch Method
    case 'lm'
      dtheta = lm_step(J,err,delta);
      pc_new = AdjustParameters(pc,dtheta);
    case 'sls'
      dtheta = lm_step(J,err,delta); % sls not yet implemented
      pc_new = AdjustParameters(pc,dtheta); % this will need LS solution too.
  end
  
  y_est = nlsim(pc_new,u);
  new_err = y - y_est;
  vaf_new = double(vaf(y,y_est));
  if vaf_new > vaf_old
    pc = pc_new;
    update_jacobian = 1;
    improve = (vaf_new - vaf_old);
    vaf_old = vaf_new;
    delta = delta*accel;
    err = new_err;
    disp([count vaf_old improve]);
    if improve < Threshold
      count = MaxIts;
    end
  else
    update_jacobian = 0;
    delta = delta*decel;
    disp([count vaf_old]);
  end
  count = count + 1;
end

set(pc,'comment',...
    'identified using Levenberg-Marquardt Optimization');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dtheta = lm_step(J,err,delta);
% Levenberg-Marquard step, implemented as in Sjoberg and Viberg (1997).


[N,NPar] = size(J);
dtheta = [J;sqrt(delta)*eye(NPar)]\[double(err);zeros(NPar,1)];


function pc = AdjustParameters(pc,dtheta);
% apply a change dtheta to the parameters in a parallel casacde model


P = get(pc,'parameters');
assign(P);
stuff = get(pc,'elements');

% update IRFs
for i = 1:NPaths
  offset = NLags*(i-1);
  h = stuff{i,1};
  hd = get(h,'data');
  hd = hd + dtheta(offset+1:offset+NLags);
  set(h,'data',hd);
  stuff{i,1} = h;
end  

% update first nonlinearity (including zero order term);
m = stuff{1,2};
coeffs = get(m,'coef');
offset = NLags*NPaths;
coeffs = coeffs + dtheta(offset+1:offset+OrderMax+1);
set(m,'coef',coeffs);
stuff{1,2} = m;

% update remaining nonlinearities (no zero order terms);
for i = 2:NPaths
  offset = NLags*NPaths + OrderMax*(i-1)+1;
  m = stuff{i,2};
  coeffs = get(m,'coef');
  coeffs = coeffs + [0;dtheta(offset+1:offset+OrderMax)];
  set(m,'coef',coeffs);
  stuff{i,2} = m;
end

set(pc,'elements',stuff);



function J = Jacobian(pc,u)
% computes Jacobian matrix for parallel cascade model.
% column order, all IRF weights followed by the zero-order polynomial
% coefficient, followed by the first and higher-order polynomial
% coefficients, path by path 


P = get(pc,'parameters');
assign(P);

NPar = NPaths*(OrderMax+NLags)+1;
N = length(u);
J = ones(N,NPar); % sets up the zero order polynomial coefficient column.


stuff = get(pc,'elements');
h = stuff{1,1};
Ts = get(h,'domainincr');


for i = 1:NPaths
  h = stuff{i,1};
  x = nlsim(h,u);
  m = stuff{i,2};
  mprime = ddx(m);
  mpx = Ts*double(nlsim(mprime,x));
  udel = double(u);
  
  offset = NLags*(i-1);
  for j = 1:NLags
    J(:,offset+j) = udel.*mpx;
    udel = [0;udel(1:N-1)];
  end
  
  offset = NLags*NPaths + 1 + (i-1)*OrderMax;
  mtest = m;
  for i = 1:OrderMax
    coeffs = zeros(OrderMax+1,1);
    coeffs(i+1) = 1;
    set(mtest,'coef',coeffs);
    J(:,offset+i) = double(nlsim(mtest,x));
  end
    
end   

    


function pc = init_random(pc,Ts);
% random initialization of pc model


P = get(pc,'parameters');
assign(P);

stuff = get(pc,'elements');
h = stuff{1,1};
m = stuff{1,2};

stuff = cell(NPaths, 2);

hdata = zeros(NLags,1);
coeffs = zeros(OrderMax+1,1);
set(h,'data',hdata,'NLags',NLags,'domainincr',Ts);
set(m,'order',OrderMax,'coef',coeffs,'type','power');

for i = 1:NPaths
  hdata = randn(NLags,1);
  set(h,'data',hdata);
  stuff{i,1} = h;
  
  coeffs(2:3) = randn(2,1);
  set(m,'coef',coeffs);
  stuff{i,2} = m;
end

set(pc,'elements',stuff);
