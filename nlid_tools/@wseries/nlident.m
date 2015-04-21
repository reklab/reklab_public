function ws = nlident (ws, z, varargin)
% Identify  a series (wiener series object)

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

if nargin > 2
   set(ws,varargin);
end

if isa(z,'pcascade'),
  ws=pc2ws (ws,z);
elseif isa(z,'lnbl')
  ws = ln2ws(z,ws);
elseif isa(z,'nlbl')
  ws = nl2ws(z,ws);
elseif isa(z,'lnlbl'),
  ws=lnl2ws(z,ws);
elseif isa(z,'vseries')
  ws = vs2ws(z,ws);
elseif isa(z,'nlm') 
  ws = nlm2ws(z,ws);
elseif isa(z,'nldat') | isa(z,'double')
  P=get(ws,'parameterSet');
  assign (P)   
  zd=double(z);
  u=zd(:,1);
  % remove the input mean.
  u = u - mean(u);
  
  % compute the input variance (used to normalize the kernels).
  uvar = std(u)^2;
  
  y=zd(:,2);
  if isnan(nLags),
    nLags=min (16,length(y)/100);
  end
  if isa(z,'nldat')
    Ts=get(z,'domainIncr');
  else
    % if the data is a double, get Ts from the wseries object.
    temp = get(ws,'elements');
    Ts = get(temp{1},'domainIncr');
  end
  
    
  % set up a cell array to accomodate the kernels.
  elements = cell(orderMax+1,1);
  switch lower(method),
    case 'ls',
      %
      % Compute kernels using Lee-Schetzen cross correlation.
      % set up an empty correlation object 
      % make it a covariance, so that input/output means are always
      % subtracted 
      phi = cor;
      set(phi,'nLags',nLags,'corType','covar');
      
      %  the zero order kernel is simply the output mean. 

      k0 = mean(y);
      wk0 = wkern;
      set(wk0,'domainIncr',Ts,'kernOrder',0,'dataSet',k0);
      elements{1} = wk0;
      
      
      if orderMax > 0
	%  the first-order kernel is computed from the first-order 
	%  cross-correlation between the input and residuals (after the
	%  zero-order output is  subtracted) -- divided by the input 
	%  variance.
	
	r0 = y  - k0;
	k1 = double(nlident(phi,[u r0],'kernOrder',1))/(Ts*uvar);
	wk1=wkern(wk0,'kernOrder',1,'dataSet',k1);
	elements{2} = wk1;
      end

      if orderMax > 1
	%  The output of the first-order kernel is subtracted, and the
	%  second-order kernel obtained from the second-order
	%  cross-correlation -- divided by the input variance squared.
	%  Also note the multiplication by a factor of 0.5.
	
	r1 = r0 - Ts*filter(k1,1,u);
	k2 = double(nlident(phi,[u r1],'kernOrder',2))/(2*uvar^2*Ts^2);
	wk2 = wkern(wk0,'kernOrder',2,'dataSet',k2);
	elements{3} = wk2;
      end
      
      if orderMax > 2
	r2 = r1 - double(nlsim(wk2,u,uvar));
	k3 = double(nlident(phi,[u r2],'kernOrder',3))/(6*uvar^3*Ts^3);
	wk3 = wkern(wk0,'kernOrder',3,'dataSet',k3);
	elements{4} = wk3;
      end
      
      if orderMax > 3
	warning('fourth (and higher) order kernels not estimated');
      end
      
      
      
    case 'toeplitz',
      numsides=1;
      [k0,k1,k2]=wienk(u,y,nLags,numsides);
  
      k1=k1/Ts;
      k2=k2/Ts^2;
      wk0=wkern('domainIncr',Ts, 'kernOrder',0,...
	  'comment','Zero order Wiener kernel', 'dataSet',k0);
      wk1=wkern('domainIncr', Ts,'kernOrder',1,...
	  'comment','First order Wiener kernel', 'dataSet',k1);
      wk2=wkern(k2,'domainIncr',Ts,'kernOrder',2,...
	  'comment','Second order Wiener kernel','dataSet',k2);
      elements = {wk0;wk1;wk2};
  end
  
  



  set(ws,'elements',elements,'variance',uvar,'nLags',nLags);

  
else
  error (['wseries cannot handle input of Class:' class(z)]);
end
%	set (ws,'Comment',[' Wiener series identified by ' Method ]);
return
% @wseries/nlident



function ws=pc2ws(wsin,pc);
ordermax=get(wsin,'orderMax');
k=cell(ordermax+1,1);
for order=0:ordermax
   wskern = pcas2wiener(pc,order);
   k{1,order+1}=wskern;
end
ws=wsin;
set(ws,'elements', { k{1}; k{2}; k{3}});


