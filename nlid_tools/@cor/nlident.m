function C = nlident (Cin, A, varargin);
% cor/nlident - identify correlation objects
%

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

C=Cin;
A=nldat(A);
if nargin < 2,
   disp('nlident takes two inputs for cor objects: cor, Z' );
elseif nargin > 2,
   set(C,varargin{:});
end
tempComment=C.comment;
[nsamp,nchan,nreal]=size(A);
incr=A.domainIncr;
domainName=A.domainName;
if isnan(incr),
   error('cor estimation requires domainIncr to be specified');
end

ninputs=1;
noutputs=nchan-ninputs;
assign(C.parameterSet);

if isnan(nLags),
   nLags = nsamp/10;
end
% Remove mean for covariance and coefficient functions
%
if any(strcmp(corType,{'covar' 'coeff'})),
   A=A-mean(A);
end

aDouble=double(A);
x=aDouble(:,1);

switch noutputs
  case 0
    % auto-correlation, so use input as output and compute cross
    % correlation.
    y = x;
  case 1
    % cross-correlation, output is second column of Adata
    y = aDouble(:,2);
  otherwise
    % we have multiple outputs, so bail
    error ('multiple outputs not yet implemented');
end

switch kernOrder
  case 0
    cx = mean(y);
    
  case 1
    [xc]=xcorr(y,x,round(nLags-1),biasMode);
    if nSides == 1
      xc = xc(nLags:end);
      start = 0;
    else
      start= -(nLags-1)*incr;
  end
  
case 2
  if nSides == 1
    yd = y;
    M = nLags;
    start = 0;
  else
    yd = zeros(nsamp,1);
    yd(nLags:end) = y(1:nsamp-nLags+1);
    start = -incr*(nLags-1);
    M = 2*nLags-1;
  end
  xc=corx2y(x,yd,M);
  biasMode = 'biased';
  
case 3
  if nSides == 1
    yd = y;
    M = nLags;
    start = 0;
  else
    yd = zeros(nsamp,1);
    yd(nLags:end) = y(1:nsamp-nLags+1);
    start = -incr*(nLags-1);
    M = 2*nLags-1;
  end
   xc=corx3y(x,yd,M);
% warning ('corrx3y not working');
 xc = reshape(xc,M,M,M);
xc=zeros(M,M,M);
  biasMode = 'biased';
    
otherwise
  error('fourth and higher order correlations not yet implemented');
end
  
if strcmp(corType,'coeff')
  % coefficient function, so normalize.
  xc = xc/(std(x)^kernOrder*std(y));
end

set(C,'dataSet',xc,...
   'nLags',nLags', ...
   'nSides',nSides,...
   'domainIncr',incr, ...
   'biasMode',biasMode,...
   'domainStart',start, 'domainName',domainName);
if strcmp('Default comment', tempComment),
    set(C,'comment',['Correlation for:' A.comment]);
else
    set(C,'comment',tempComment);
end

end

% ... cor/nlident  

