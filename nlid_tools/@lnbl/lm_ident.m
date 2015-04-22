function bl = lm_ident(bl,iodata)
% nonlinear optimization method for Wiener cascade fitting.
% uses a Levenberg-Marquardt second-order gradient descent search.
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 



subsys = get(bl,'elements');
h_est = subsys{1};
m_est = subsys{2};

u = iodata(:,1);
z = iodata(:,2);

MaxIts = get(bl,'nMaxIts');
Thresh = get(bl,'threshNMSE');
Accel = get(bl,'accel');
Decel = get(bl,'decel');
Delta = get(bl,'delta');


y_est = nlsim(bl,u);
err = z - y_est;
vaf_old = double(vaf(z,y_est));

count = 1;
update_jacobian = 1;
bl_new = bl;

while count <= MaxIts
  if update_jacobian
    J = wcas_jacobian(h_est,m_est,u);
  end
  [h_new,m_new] = lm_step(h_est,m_est,J, err, Delta);
  set(bl_new,'elements',{h_new,m_new});
  y_est = nlsim(bl_new,u);
  new_err = z - y_est;
  vaf_new = double(vaf(z,y_est));
  if vaf_new > vaf_old
    h_est = h_new;
    m_est = m_new;
    update_jacobian = 1;
    improve = (vaf_new - vaf_old)/100;
    vaf_old = vaf_new;
    Delta = Delta*Accel;
    err = new_err;
    if improve < Thresh
      count = MaxIts;
    end
  else
    update_jacobian = 0;
    Delta = Delta*Decel;
  end
  count = count + 1;
end

set(bl,'elements',{h_est,m_est},'comment',...
    'LN model identified using Levenberg-Marquardt Optimization');



function [hnew,mnew] = lm_step(h_old,m_old, J, err, delta);

h = get(h_old,'dataSet');
m = get(m_old,'polyCoef');

theta = [h;m];
M = length(theta);

% check this implementation vs. Sjoberg and Viberg
Je = [J;sqrt(delta)*eye(M)];
erre = [double(err);zeros(M,1)];
theta_new = theta + Je\erre;

%H = (J'*J + delta*eye(M));
%theta_new = theta + inv(H)*J'*err;

hnew = h_old;
hlen = get(hnew,'nLags');
set(hnew,'dataSet',theta_new(1:hlen));
mnew = m_old;
set(mnew,'polyCoef',theta_new(hlen+1:end));



function J = wcas_jacobian(h,m,u);

hh = get(h,'dataSet');
Ts = get(h,'domainIncr');
hlen = length(hh);

Q = get(m,'polyOrder');



x = nlsim(h,u);
mprime = ddx(m);
mpx = Ts*double(nlsim(mprime,x));


udel = double(u);
N = length(udel);
J = zeros(N,hlen+Q+1);


for i = 1:hlen
  J(:,i) = udel.*mpx;
  udel = [0;udel(1:N-1)];
end

mtest = m;
for i = 0:Q
  coeffs = zeros(Q+1,1);
  coeffs(i+1) = 1;
  set(mtest,'polyCoef',coeffs);
  J(:,i+hlen +1) = double(nlsim(mtest,x));
end




  
  
  


