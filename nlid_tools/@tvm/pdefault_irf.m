function IOut = pdefault(IIn);
% set default parameters for irfs

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

%p=get(IIn,'parameters');
%j=size(p);
j=1;
p=param('name','NLags','default',16,'help','Number of lags', ...
   'type','real','Limits', [ 1 inf]);
p(2)=param('name','NSides','default',1,'help','number of sides', ...
   'type','real','limits',[1 2]);
p(3)=param('name','Order','default',1,'help','kernel order', ...
   'type','real','limits',[1 2]);
p(4)=param('name','Level','default',95,'help','error level');
p(5)=param('name','Mode','default','full','help','mode ');
p(6)=param('name','Method','default','tvfil','help','Identification method: tvfil/corr/pseduo');
p(7)=param('name','DisplayFlag','default','true','help','display');
IOut=IIn;
set(IOut,'Parameters',p);
