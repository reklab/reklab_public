function IOut = pdefault(IIn);
% set default parameters for vseries

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

p=get(IIn,'parameters');
j=size(p);

method = get(IIn,'Method');
if j>1
  nlags = get(IIn,'NLags');
else
  nlags = NaN;
end

% method is initially set to FOA in vseries/mkvseries.   The next line 
% sets p(1) to contain method, and it defaults to whatever it used to be.
% the idea is to estabish the help, type and limits properties.
%
% If p has more than 1 element, NLags has been defined.  Therefore, it
% should not be changed, since this parameter is valid regardless of the
% identification method.  Thus, if it has not been previously set, NLags
% defaults to NaN, however, once set, it remains sticky. 

p(1) = param('name','Method','default',method,...
    'help','Identification method' , 'type','select',...
    'limits', {'FOA' 'foa' 'Laguerre' 'laguerre'});
p(2)=param('name','NLags','default',nlags,...
    'help','Number of lags in each kernel' ,...
   'type','number','limits', {0 1000});
p(3)=param('name','OrderMax','default',2,'help','Maximum order for series' ,...
   'type','number','limits', {0 3});

switch lower(method)
  case 'foa'
    % set limits on OrderMax to [0 2]
    % since FOA is limited to second-order systems
    p = p(1:3);
    p(3)=param('name','OrderMax','default',2,...
	'help','Maximum order for series' ,...
	'type','number','limits', {0 2});
  case 'laguerre'
    p(4) = param('name','alpha','default',NaN,...
	'help','decay parameter for Laguerre filters',...
	'type','number','limits',{eps 1-eps});
    p(5) = param('name','NumFilt','default',NaN,...
	'help','Number of Laguerre filters in expansion',...
	'type','number','limits',{1 100});
    if isnan(nlags)
      max_delay = 1000;
    else
      max_delay = nlags-1;
    end
    p(6) = param('name','delay','default',0,...
	'help','delay before start of Laguerre fitlers',...
	'type','number','limits',{0 max_delay});
    p=p(1:6);
  otherwise
    error('unrecognized identification method');
end

IOut=IIn;
set(IOut,'Parameters',p);
