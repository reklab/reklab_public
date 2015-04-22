function ws=nl2ws (nl, wsin)
% generate Wiener kernels for Wiener cascade model

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


% get the input variance from the provided wseries model
sigma2 = get(wsin,'variance');

stuff = get(nl,'elements');
% first extract the linear subsystem
h = double(stuff{2});
Ts = get(stuff{2},'domainIncr');
hlen = length(h);

% now deal with the nonlinearity
p = stuff{1};
if isnan(sigma2) 
  % get the variance from the polynomial
  ptemp = nlident(p,'polyType','hermite');
  sigma2 = get(ptemp,'polyStd')^2;
  if isnan(sigma2)
    sigma2 = 1;
  end  
end



% convert the nonlinearity into a hermite polynomial, normalized for an
% input with Std sigma_x
polytype = get(p,'polyType');
if strcmp(polytype,'hermite')
  p = nlident(p,'type','tcheb');
end
sigma = sqrt(sigma2);
set(p,'polyStd',sigma,'polyMean',0);

ph = nlident(p,'polyType','hermite');
coeffs = get(ph,'polyCoef');

orderMax=get(wsin,'orderMax'); 
Q = min(length(coeffs), orderMax+1);
kernels = cell(Q,1);


for q = 1:Q,
  den = Ts^(q-2)*sigma^(q-1);
  k = (coeffs(q)/den)*hckern(h,q-1);
  kernels{q} = wkern('kernOrder',q-1,'domainIncr',Ts, 'dataSet',k);
end
ws=wsin;
set(ws,'elements',kernels,'variance',sigma2,'orderMax',Q-1,'nLags',hlen,...
    'comment','Transformed Hammerstein Cascade');
