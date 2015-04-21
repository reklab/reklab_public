function pc = pcm_nlls(pc,z);
% nonlinear least squares optimization methods to fit parallel cascade
% models.  Normally called from nlident.m


% 1 extract data and tuning parameters.

u = z(:,1);
y = z(:,2);
N = length(y);
Ts = get(u,'domainIncr');

P = getParamValStruct(pc.parameterSet);

% creates variables containing parameters:
% Method, NPaths, NLags, OrderMax, MaxIts, Threshold, 
% accel, decel, delta, initialization  

NPaths=P.nPaths;
if isnan(NPaths)
% default setting is NPaths = NaN.  This works for parallel cascade
% algorithms, as NPaths can be determined on the fly. 
% Since NLLS algorithms require a fixed structure, we will have to set
% a default.  This computation gives at least 10 data points per parameter.
  NPaths = floor(N/(10*(P.nLags+P.polyOrderMax)));
end

% initialize with random IRF weights and polynomial coefficients, if
% requested. 
if strcmp(P.initMethod,'random')
  pc = init_random(pc,Ts);
end


y_est = nlsim(pc,u);
err = y - y_est;
vaf_old = double(vaf(y,y_est));

count = 1;
update_jacobian = 1;
pc_new = pc;
 delta=P.searchDelta;
while count <= P.nMaxIts
  if update_jacobian
    J = Jacobian(pc,u);
  end
  switch P.idMethod
    case 'lm'
      dtheta = lm_step(J,err,P.searchDelta);
      pc_new = AdjustParameters(pc,dtheta);
    case 'sls'
      disp('sls not implemented so use lm');
      dtheta = lm_step(J,err,P.searchDelta); % sls not yet implemented
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
    delta = delta*P.searchAccel;
    err = new_err;
    disp(['Count:' num2str(count) ' VAF:' num2str(vaf_old) ' Improvement:' num2str(improve) ]);
    if improve < P.searchThreshold
      count = P.nMaxIts;
    end
  else
    update_jacobian = 0;
    delta = delta*P.searchDecel;
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


P =getParamValStruct(pc.parameterSet); 
stuff = pc.elements;

% update IRFs
for i = 1:P.nPaths
  offset = P.nLags*(i-1);
  h = stuff{i,1};
  hd = get(h,'dataSet');
  hd = hd + dtheta(offset+1:offset+P.nLags);
  set(h,'dataSet',hd);
  stuff{i,1} = h;
end  

% update first nonlinearity (including zero order term);
m = stuff{1,2};
coeffs = get(m,'polyCoef');
offset = P.nLags*P.nPaths;
coeffs = coeffs + dtheta(offset+1:offset+P.polyOrderMax+1);
set(m,'polyCoef',coeffs);
stuff{1,2} = m;

% update remaining nonlinearities (no zero order terms);
for i = 2:P.nPaths
  offset = P.nLags*P.nPaths + P.polyOrderMax*(i-1)+1;
  m = stuff{i,2};
  coeffs = get(m,'polyCoef');
  coeffs = coeffs + [0;dtheta(offset+1:offset+P.polyOrderMax)];
  set(m,'polyCoef',coeffs);
  stuff{i,2} = m;
end

set(pc,'elements',stuff);



function J = Jacobian(pc,u)
% computes Jacobian matrix for parallel cascade model.
% column order, all IRF weights followed by the zero-order polynomial
% coefficient, followed by the first and higher-order polynomial
% coefficients, path by path 


P = getParamValStruct(pc.parameterSet); 


NPar = P.nPaths*(P.polyOrderMax+P.nLags)+1;
N = length(u);
J = ones(N,NPar); % sets up the zero order polynomial coefficient column.


stuff = pc.elements;
h = stuff{1,1};
Ts = h.domainIncr;


for i = 1:P.nPaths
  h = stuff{i,1};
  x = nlsim(h,u);
  m = stuff{i,2};
  mprime = ddx(m);
  mpx = Ts*double(nlsim(mprime,x));
  udel = double(u);
  
  offset = P.nLags*(i-1);
  for j = 1:P.nLags
    J(:,offset+j) = udel.*mpx;
    udel = [0;udel(1:N-1)];
  end
  
  offset = P.nLags*P.nPaths + 1 + (i-1)*P.polyOrderMax;
  mtest = m;
  for i = 1:P.polyOrderMax
    coeffs = zeros(P.polyOrderMax+1,1);
    coeffs(i+1) = 1;
    set(mtest,'polyCoef',coeffs);
    J(:,offset+i) = double(nlsim(mtest,x));
  end
    
end   

    


function pc = init_random(pc,Ts);
% random initialization of pc model


P = getParamValStruct(pc.parameterSet);


stuff = pc.elements;
h = stuff{1,1};
m = stuff{1,2};

stuff = cell(P.nPaths, 2);

hdata = zeros(P.nLags,1);
coeffs = zeros(P.polyOrderMax+1,1);
set(h,'dataSet',hdata,'nLags',P.nLags,'domainIncr',Ts);
set(m,'polyOrder',P.polyOrderMax,'polyCoef',coeffs,'polyType','power');

for i = 1:P.nPaths
  hdata = randn(P.nLags,1);
  set(h,'dataSet',hdata);
  stuff{i,1} = h;
  coeffs(2:3) = randn(2,1);
  set(m,'polyCoef',coeffs);
  stuff{i,2} = m;
end

set(pc,'elements',stuff);
