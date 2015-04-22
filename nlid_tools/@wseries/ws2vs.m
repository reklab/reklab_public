function vs = ws2vs(ws,vsin);
% converts a wiener series into a volterra series
% uses theorem 5.3 from Rugh (1981).

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


wkernels = get(ws,'elements');
Ts = get(wkernels{1},'domainIncr');
sigma2 = get(ws,'variance');
Nlags = get(ws,'nLags');
Q = get(ws,'orderMax');

vkernels = cell(Q+1,1);

for q = Q:-1:0   
    k_volt = vkern('domainIncr',Ts,'kernOrder',q, 'dataSet',zerokern(Nlags,q));
  for i = 0:2:Q-q
    j = i/2;
    num = ((-1)^j)*factorial(q+i)*sigma2^j;
    den = factorial(q)*factorial(j)*2^j;
    gain = num/den;
    kred = wkernels{q+i+1};
    for l = 1:j
      kred = reduce_kernel(kred)*Ts^2;
    end
   
    newKern= double(k_volt) + gain*double(kred);
    set (k_volt,'dataSet',newKern); 
  end
  vkernels{q+1}  = k_volt;
end

vs = vsin;
set(vs,'elements',vkernels);


%% a much more efficient solution would be to generate all of the reduced
%% kernels, and store them in a cell array, as this would avoid regenerating
%% the same reduced kernels, over and over again. 


function k = reduce_kernel(kern)
% integrates kern over the last pair of dimensions

old_order = get(kern,'kernOrder');
new_order = old_order - 2;
Nlags = get(kern,'nLags');

k = vkern;
set(k,'kernOrder',new_order,'nLags',Nlags);


kdat = zerokern(Nlags,new_order);
kk = double(kern);
switch old_order
  case 0 
    error('oops, knocking 2 dimensions off of an order 0 kernel');
  case 1 
    error('oops, knocking 2 dimensions off of a first order kernel');
  case 2
    kdat = sum(diag(kk)); 
  otherwise  
    command = 'kdat = kdat + squeeze(kk(:';
    for n = 2:old_order-2
      command = [command ',:'];
    end
    command = [command ',m,m));'];
      for m = 1:Nlags
	eval(command);
      end
  end
set(k,'dataSet',kdat);  


