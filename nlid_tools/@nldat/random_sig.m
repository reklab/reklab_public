function Y = random_sig (X , C, NCHANGE, NREPORT );
%
% X     - seed signal 
% C     - desired autocorelation function
% NCHANGE= number of interchanges
% NREPORT -- frequency to report
%

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if (nargin <4),
 NREPORT=NCHANGE/10;
end
x=X.Data;
cdom=domain(C);
izero=find(cdom==0);
cdata=get(C,'data');
c=cdata(izero:end)';
y=random_sig_exec( x,c, NCHANGE, NREPORT);
Y=nldat(y,'domainincr',X.DomainIncr);
return;




% 
% Random sig wrapper for NLID toolbox



function [y,cout] = random_sig_exec ( x, c, nchange, nreport );
% MATLAB implementation of the NEXUS function SIG;
% Generate noise with specified  PDF and spectra
%          x       - seed signal
%          c       - desired covariance function for lags >= 0; 
%          nchange - number of interchanges
%          nreport - report frequency
% 5 Apr 94 rek
% 2002 Oct 11 Normalize input covariance function
if (nargin <4),
 nreport=nchange/10;
end
icount =0;
n = length(x);
kmax = length(c);
y=x-mean(x);
% normalize input covariance
c=c/c(1); 
cov1=double(cor(y,'NLags',kmax,'bias','biased','type','covar'))';
c1 = cov1(1);
cov2(1)=c1;
ss1 = cov(c-(cov1/c1));
ss_ref=ss1;
subplot (1,1,1);
plot (0,ss1/ss_ref,'o');
axis ([ 0 nchange 0 1]);
drawnow
hold
for i=1:nchange,
  % {{{ Generate a possible interchange

  j= ceil(rand (2,1)*n);
  i1 = j(1);
  i2 = j(2);
  x1 = y(i1);
  x2 = y(i2);

% }}}
  % {{{ Update  autocorrelation estimate

 
  for k=2:kmax,
    sum1=0;
    sum2=0;
    if (i1+k-1 == i2),
      if(i1-k+1 >= 1),
	sum1 = sum1 + x1 * y(i1-k+1);
	sum2 = sum2 + x2 * y(i1-k+1);
      end
      if (i2+k-1 <= n),
	sum1 = sum1 + x2 * y(i2+k-1);
	sum2 = sum2 + x1 * y(i2+k-1);
      end
    elseif (i1-k+1 == i2),
      if (i1+k-1 == i2),
	sum1 = sum1 + x1 * y(i1+k-1);
	sum2 = sum2 + x2 * y(i1+k-1);
      end
      if (i2-k+1 >=1),
	sum1 = sum1 + x2 * y(i2-k+1);
	sum2 = sum2 + x1 * y(i2-k+1);
      end
    else,
      if (i1+k-1 <=n),
	sum1 = sum1 + x1 * y(i1+k-1);
	sum2 = sum2 + x2 * y(i1+k-1);
      end
      if (i1-k+1 >=1),
	sum1 = sum1 + x1 * y(i1-k+1);
	sum2 = sum2 + x2 * y(i1-k+1);
      end
      if (i2+k-1 <=n),
	sum1 = sum1 + x2 * y(i2+k-1);
	sum2 = sum2 + x1 * y(i2+k-1);
      end
      if (i2-k+1 >=1),
	sum1 = sum1 + x2 * y(i2-k+1);
	sum2 = sum2 + x1 * y(i2-k+1);
      end
    end
    cov2(k) = (cov1(k)*n +sum2 -sum1)/n ;
  end

% }}}
  % {{{ Check for improvement

ss2 = cov(c-(cov2/c1));
if (ss2 < ss1),
  y(i1)=x2;
  y(i2)=x1;
  cov1=cov2;
  ss1=ss2;
  end

% }}}
  % {{{ Monitor progress

  if (i > icount),
    icount = icount+nreport;
    disp([ int2str(i) ' ' num2str(ss1/ss_ref) ]);
    subplot (2,1,1);
    plot (i,ss1/ss_ref,'o');
    subplot (2,1,2);
    plot ([c; cov2/c1]');
    drawnow;
  end

% }}}
end
hold off
cout = cov1/c1;
% {{{ Emacs local variables

% Local variables:
% folded-file: t
% end;

% }}}







