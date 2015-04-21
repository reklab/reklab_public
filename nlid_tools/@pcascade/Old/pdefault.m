function IOut = pdefault(IIn);
% set default parameters for pcascade

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

p=get(IIn,'parameters');
method = get(IIn,'Method');
% Method is initially set in pcascade/mkpc to be 'eig'.  The first line
% sets p(1) to contain method, and it defaults to whatever it used to be.
% the idea is to estabish the help, type and limits properties.
j=size(p);
p(1)=param('name','Method','default',method,'help','Identification method' ,...
    'type','select','limits', {'eig' 'gen_eig' 'lm' 'slice' 'sls'});
if j < 2
  p(2) = param('name','NPaths','default',NaN,...
      'help','Number of paths','type','number',...
      'limits',{1 1000});
end
if j<3
  p(3)=param('name','NLags','default',16,'help',...
      'Number of lags in each kernel' ,'type','number','limits', {0 1000}); 
end
if j <4
  p(4)=param('name','OrderMax','default',5,'help',...
      'Maximum order for series' ,'type','number','limits', {0 10});
end

switch lower(method)
  case{'eig', 'gen_eig','slice'}
    % if previous method was eig, gen_eig or slice, then mode has already
    % been set, so preserve its current value
    set_mode = 0;
    if j < 5
      set_mode = 1;
    elseif ~strcmp(get(p(5),'name'),'mode')
      set_mode = 1;
    end
    if set_mode
    p(5) = param('name','mode','default','auto','help',...
	'method to select polynomial order','type','select',...
	'limits',{ 'auto' 'full' 'manual' });
    end
end    


switch lower(method)
  case {'eig', 'gen_eig'}
    p = p(1:5);
  case 'slice'
    p(6) = param('name','MaxPaths','default',10,...
	'help','Maximum number of paths tested in search',...
	'type','number','limits',{1 100});
    p(7) = param('name','TestPaths','default','yes',...
	'help','Retain only statistically significant paths',...
	'type','select','limits',{'yes','y','no','n'});
    p = p(1:7);
  case {'lm', 'sls'}
    p(5) = param('name','MaxIts','default',100,...
	'help','Maximum number of paths tested in the expansion',...
	'type','number','limits',{1 10000});
    p(6)=param('name','Threshold','default',.01,'help',...
	'NMSE for success','type','number','limits',{0 1});
    p(7)=param('name','accel','default',.8,'help',...
	'ridge multiplied by accell after successful update',...
	'type','number','limits', {0.001 0.999});
    p(8)=param('name','decel','default',2,'help',...
	'ridge multipled by devel after unsuccessful update',...
	'type','number','limits', {1.0001 inf});
   p(9)=param('name','delta','default',10,'help',...
       'initial size of ridge added to Hessian','type','number',...
       'limits',{0 inf});
    p(10)= param('name','initialization','default','random',...
	'help','Starting model specification','type','select',...
	'limits',{'current','random'});
  otherwise    
    error(['pcascade objects do not support method: ' method]);
end




IOut=IIn;
set(IOut,'Parameters',p);
