function vs = lag_id(vs,u,y,Ts)
%  Laguerre based identification

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


Qmax = get(vs,'vsOrderMax');
hlen = get(vs,'nLags');
alpha = get(vs,'alphaLaguerre');
num_filt = get(vs,'nLaguerreFilt');
delay = get(vs,'delayLaguerre');


kernels = get(vs,'elements');
vk0 = kernels{1};
set(vk0,'domainIncr',Ts,'nLags',hlen);

% set any missing parameters
if isnan(num_filt)
  if isnan(alpha)
    num_filt = ceil(hlen/4);
  else
    num_filt = num_filters(hlen,alpha,delay);
  end
  set(vs,'nLaguerreFilt',num_filt);
end

if isnan(alpha)
  alpha = ch_alpha(hlen,num_filt,delay);
  set(vs,'alphaLaguerre',alpha);
end

  
u_basis = laguerre_basis(u,num_filt,alpha,delay);

static = polynom;
set(static,'polyOrderMax',Qmax,'nInputs',num_filt,'polyType','power','polyOrderSelectMode','full');
static = nlident(static,[u_basis y]);

coeff = get(static,'polyCoef');


impulse = [1;zeros(hlen-1,1)];
basis = laguerre_basis(impulse,num_filt,alpha,delay);



kernels = cell(Qmax+1,1);

h0 = coeff(1);
kernels{1} = vkern(vk0,'kernOrder',0,'dataSet',h0);



for i = 1:Qmax
  ki = gen_kern(basis,coeff,i);
  kernels{i+1} = vkern(vk0,'kernOrder',i,'dataSet',ki/(Ts^i));
end

set(vs,'elements',kernels);
set(vs,'comment','Identified using Laguerre Expansion Technique');

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


function nfilt = num_filters(hlen,alpha,delay)

imp = [1;zeros(2*hlen,1)];
basis = laguerre_basis(imp,hlen,alpha,delay);

searching = 1;
nfilt = 1;

while searching
  idx = max(find(abs(basis(:,nfilt))>0.01));
  if idx > hlen-1
    searching = 0;
    nfilt = nfilt -1;
  else
    nfilt = nfilt + 1;
  end
end

%test = abs(basis(end,:));
%nfilt = max(find(test<0.01));


function alpha = ch_alpha(hlen,num_filt,delay);

alpha = 0.5;
ntest = num_filters(hlen,alpha,delay);
step = 0.25;
if ntest == num_filt
  searching = 0;
else
  searching = 1;
end

while searching
  if ntest < num_filt
    alpha = alpha - step;
  else
    alpha = alpha + step;
  end
  ntest = num_filters(hlen,alpha,delay);
  if ntest == num_filt
    searching = 0;
  else
    step = step/2;
  end
end

