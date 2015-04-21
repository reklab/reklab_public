function f = random_pdf (type, x, p1, p2,p3)
% generates proability densities for various distirbutions
%
% usage:
%          f = random_pdf (type, x, p1, p2, p3 )
%
% where
%     x are points to evaluate pdf
%     the other parameters are used as follows:
%                  type              p1       p2         p3     
%               'normal',            mean     standard-deviation 
%               'chi-square',        dof
%               'student-t',         dof
%               'F-distribution'     dof1     dof2
%               'log-normal',        mean     variance   offset  
%               'rectangular',       xmin     xmax
%               'exponential',       mean
%               'binomial',          n        probabilty
%               'Poisson'            mean 
%

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% {{{ Normal
if (strcmp(type,'normal')),
  %
  % Parse input and check its validity
  %
  if (nargin < 4),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''normal'',x,mean,variance');
    f=NaN;
    return
  end
  mean =p1;
  sd = p2;
  if (sd <= 0),
    disp('ERROR: standard deviation must be >0');
    f=NaN;
    return
  end
    
  %
  % generate pdf
  %
  x1 = (x-mean).^2;
  x2 = x1./(2*sd*sd);
  f = exp(-x2)./(sd*sqrt(2*pi));

% }}}
% {{{ chi-square

elseif(strcmp(type,'chi-square')),
  %
  % parse input and check its validity
  %
  if (nargin < 3),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''chi-square'',x,degrees_of_freedom');
    f=NaN;
    return
  end
  if (p1 <0),
    disp('ERROR: degrees-of-freedom must be > 0 for chi-square');
      f=NaN;
      return
    end
    if (any(x<0)),
    disp('ERROR: x >0 for chi-square')
    f=NaN;
    return
  end
 %
 % Gnerate pdf
 % 
 c1=((2^(p1/2)) * gamma(p1/2)) ; 
  f1 = exp(-x/2) ;
  f2 = x.^((p1/2)-1);
  f = (f1.*f2)/c1;

% }}}
% {{{ Student-t

elseif(strcmp(type,'student-t')),
  % {{{  Parse input and check its validity

  %
  if (nargin < 3),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''student-t'',x,degrees_of_freedom');
    f=NaN;
    return
  end
  dof =p1;
  if (dof <= 0),
    disp('ERROR: degrees-of-freedom  must be >0');
    f=NaN;
    return
  end

% }}}
  % {{{ Generate pdf

  c1 = 1/sqrt(dof*pi);
  c2=gammaln((dof+1)/2) - (gammaln(dof/2));
  
  f1 = (1+((x.*x)/dof)).^(-(dof+1)/2);
  
  f = c1*exp(c2)*f1;
  

% }}}

% }}}
% {{{ F-distribution

elseif(strcmp(type,'F-distribution')),
  % {{{  parse input and check validity

  %
  if (nargin < 4),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''F-distirbution'',x,n1,n2');
    f=NaN;
    return
  end
  if (p1<=0) | (p2<=0),
    disp('ERROR: Both values for degrees-of-freedom must be > 0');
    f=NaN;
    return
    end
    
  if any(x<=0),
    disp('ERROR: random_pdf: Bad value. x must be > 0')
    f = NaN;
    return
  end 

% }}}
  % {{{ Generate pdf

  c1 = gammaln((p1+p2)/2) - (gammaln(p1/2) + gammaln(p2/2));
  c1 = exp(c1);
  c2 = (p1/p2)^(p1/2);
  c3=c1*c2;
  f1 = x.^((p1-2)/2);
  f2 = (1 + (p1/p2)*x).^((p1+p2)/2);
  f = c3.*(f1./f2);

% }}}

% }}}
% {{{ log-normal distribution

elseif (strcmp(type,'log-normal')),
  
  mean = p1
  var = p2
  offset=p3
  if (nargin < 5),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''log-normal'',x,mean,var,offset');
    f=NaN;
    return
  end
  
  if any(x <= offset),
    disp('ERROR: random_pdf: Bad value. x must be > offset')
    f = NaN;
    return
  end
  x1=x-offset;
  f1 = 1./(x1.*sqrt(var)*sqrt(2*pi));
  f2 = ((((log(x1)-mean).^2)./(var)));
  f3 = exp(-f2/2);
  f=f1.*f3;
  

  % }}}
% {{{ Rectangular

elseif (strcmp(type,'rectangular')),
  %
  % Parse input and check its validity
  %
  if (nargin < 4),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''rect'',x,xmin,xmax');
    f=NaN;
    return
  end
  xmin=p1;
  xmax=p2;
  %
  % generate pdf
  %
  f=zeros(size(x));
  i=find ((x<=xmax) & (x>=xmin));
  f(i)=ones(size(f(i)))./(xmax-xmin);

% }}}
% {{{ Exponential

elseif(strcmp(type,'exponential')),
  %
  % parse input and check its validity
  %
  if (nargin < 3),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''exponential'',x,mean');
    f=NaN;
    return
  end

 %
 % Gnerate pdf
 % 
 i=find (x>=0);
 j=find (x<0);
 f(j)=0;
 mu=p1;
  f(i) = (exp(-x(i)/mu))/mu;

% }}}
% {{{ Binomial

elseif(strcmp(type,'binomial')),
  % {{{  parse input and check its validity

  %
  if (nargin <4),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''binomial'',x,number-of-event,probability');
    f=NaN;
    return
  end

    if (any(x<0)),
    disp('ERROR: x >0 for binomial distribution')
    f=NaN;
    return
  end

% }}}
  % 
  % Generate pdf
  % 
  n=round(p1);
   
    prob=p2;
    x=round(x);
    f=zeros(size(x));
    index=(find(x>=0 & x <= n));
    x1=x(index);
    
    c1 = gammaln(n+1);
    c2 = gammaln(x1+1) + gammaln(n-x1+1);
    c3 = c1 - c2;
    c4 = prob .^ x1;
    c5 = (1-prob) .^ ( n - x1);
    fx = exp(c3) .* c4 .* c5;
    
    f(index)=fx;
  

  % }}}
% {{{ Poisson

elseif(strcmp(type,'Poisson')),
  %
  % parse input and check its validity
  %
  if (nargin < 3),
    disp('ERROR: random_pdf  not enough parameters specified');
    disp('Usage is: f = random_pdf(''Poisson'',x,mean');
    f=NaN;
    return
  end

    if (any(x<0)),
    disp('ERROR: x >0 for Possion distribution')
    f=NaN;
    return
  end
 %
 % Generate pdf
 % 
  lambda = p1;
  x=round(x);
  c1 = exp(-lambda);
  c2 = lambda .^ x;
  c3 = gamma(x+1);
  f = (c1 .* c2) ./ c3; 
  

  % }}}

else,
  disp('not a valid distribution');
  f=NaN;
end


% {{{ Emacs local variables

% Local variables:
% folded-file: t
% end;

    % }}}



