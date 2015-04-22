function ws=ln2ws (ln, wsin)
% generate Wiener kernels for Wiener cascade model

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


OrderMax = get(wsin,'orderMax');


% get the input variance from the provided wseries model
sigma2 = get(wsin,'variance');


stuff = get(ln,'elements');

% first extract the linear subsystem
h = double(stuff{1});
Ts = get(stuff{1},'domainIncr');
hlen = length(h);
% compute the gain applied by the linear system to the std of the input.
gain = Ts*sqrt(sum(h.^2));


% now deal with the nonlinearity
p = stuff{2};
if isnan(sigma2) 
  % get the variance from the polynomial
  ptemp = nlident(p,'polyType','hermite');
  sigma2 = (get(p,'polyStd')/gain)^2;
  if isnan(sigma2)
    sigma2 = 1;
  end  
end

% compute the Std of the input to the polynomial
sigma_x = gain*sqrt(sigma2);


% convert the nonlinearity into a hermite polynomial, normalized for an
% input with Std sigma_x
polytype = get(p,'polyType');
if strcmp(polytype,'hermite')
  p = nlident(p,'polyType','tcheb');
end
set(p,'polyStd',sigma_x,'polyMean',0);

ph = nlident(p,'polyType','hermite');
coeffs = get(ph,'polyCoef');

h_norm = h/sigma_x;

orderMax=get(wsin,'orderMax');
Q=length(ph.polyCoef); 
if Q > OrderMax
  % truncate polynomial to aviod producing high-order kernels.
  str = ['truncating Wiener series at order ' num2str(OrderMax)];
  warning(str);
  Q = OrderMax;
end
qMax=Q+1
kernels = cell(qMax,1);



for q = 1:qMax
  k = coeffs(q)*wckern(h_norm,q-1);
  kernels{q} = wkern(k,'kernOrder',q-1,'domainIncr',Ts);
end
ws=wsin;
set(ws,'elements',kernels,'variance',sigma2,'orderMax',Q,'nLags',hlen,...
    'comment','Transformed Wiener Cascade');
