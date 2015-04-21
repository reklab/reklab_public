function xc = xcorr ( A, varargin);
% overlaid xcorr for nldt objects

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

options = { {'ctype' 'coeff' 'corel\covar\coeff'} ...
    { 'bias' 'unbiased' 'Biased or unbiaised estimate'} ...
    {'maxlag' NaN 'maximum lag' } ...
    { 'Ninputs' 1 'Number of inputs' }};
if  arg_parse (options, varargin),
  return
end

incr=A.DomainIncr;
[nsamp,nchan,nreal]=size(A);
ninputs=min(ninputs,nchan);
noutputs=nchan-ninputs;
if isnan(maxlag),
  d=domain(A);
  maxlag=(max(d)-min(d))/10;
end
nmax = round (maxlag/A.DomainIncr);
if any(strcmp(ctype,{'covar' 'coeff'})),
  A=A-mean(A);
end
if noutputs==0,
  x=double(A);
  [xc,xcl]=xcorr(x,x,nmax,bias);  
  if strcmp (ctype,'coeff'),
    xc=xc/cov(x);
  end
  xc=nldat(xc);
  lagmin=min(xcl)*incr;
  set(xc,'DomainStart',lagmin,'DomainIncr',incr);
  set (xc,'Comment',['Auto:' ctype]);
elseif ninputs ==1 & noutputs ==1,
  a=double(A);
  x=a(:,1);
  y=a(:,2);
  [xc,xcl]=xcorr(x,y,nmax,ctype);  
  if strcmp(ctype,'coef'),
    xc=xc/sqrt(cov(x)*cov(y));
  end
  xc=nldat(xc);
  lagmin=min(xcl)*incr;
  set(xc,'DomainStart',lagmin,'DomainIncr',incr);
  set (xc,'Comment', ['Cross:' ctype]);
elseif npoutputs==1,
else
  error ('mupltiple outputs not yet implement');
end

% ... nldat/xcorr   
