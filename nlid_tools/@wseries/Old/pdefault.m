function IOut = pdefault(IIn);
% set default parameters for wkern

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

p=get(IIn,'parameters');
j=size(p);
method = get(IIn,'Method');
p(1) = param('name','Method','default',method,...
    'help','Identifcation Method','type','select',...
    'limits',{'ls' 'LS' 'Toeplitz' 'toeplitz'});

switch lower(method)
  case 'ls'
  p(2)=param('name','OrderMax','default',2,...
      'help','Maximum order for series' ,...
      'type','number','limits', {0 3});
case 'toeplitz'  
  p(2)=param('name','OrderMax','default',2,...
      'help','Maximum order for series' ,...
      'type','number','limits', {0 2});
end

if j < 3
  p(3)=param('name','NLags','default',NaN,...
      'help','Number of lags in each kernel' ,...
      'type','number','limits', {0 1000});
end

if j < 4
  p(4)=param('name','Variance','default',NaN,...
      'help','Input variance used in orthogonalization' ,...
      'type','number','limits', {0 inf});
end  

IOut=IIn;
set(IOut,'Parameters',p);
