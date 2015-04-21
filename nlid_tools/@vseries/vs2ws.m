function ws = vs2ws(vs,wsin);
% converts a wiener series into a volterra series
% uses theorem 5.4 from Rugh (1981).

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


vkernels = get(vs,'elements');
Ts = get(vkernels{1},'domainincr');
sigma2 = get(wsin,'variance');
if isnan(sigma2)
  sigma2 = 1;
end
Nlags = get(vs,'NLags');
Q = get(vs,'OrderMax');

wkernels = cell(Q+1,1);

for q = Q:-1:0
  k_wiener = wkern(zerokern(Nlags,q),'domainincr',Ts,'order',q);
  for i = 0:2:Q-q
    j = i/2;
    num = factorial(q+i)*sigma2^j;
    den = factorial(q)*factorial(j)*2^j;
    gain = num/den;
    kred = vkernels{q+i+1};
    for l = 1:j
      kred = reduce_kernel(kred)*Ts^2;
    end
    k_wiener = k_wiener + gain*kred;
  end
  wkernels{q+1}  = k_wiener;
end

ws = wsin;
set(ws,'elements',wkernels,'NLags',Nlags,'OrderMax',Q,'Variance',sigma2);


%% a much more efficient solution would be to generate all of the reduced
%% kernels, and store them in a cell array, as this would avoid regenerating
%% the same reduced kernels, over and over again. 


function k = reduce_kernel(kern)
% integrates kern over the last pair of dimensions

old_order = get(kern,'order');
new_order = old_order - 2;
Nlags = get(kern,'NLags');

k = wkern;
set(k,'order',new_order,'NLags',Nlags);


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
set(k,'data',kdat);  


