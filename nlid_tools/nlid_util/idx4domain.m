function idx = idx4domain (domainStart, domainIncr, domainVal )
% idx = idx4domain (domainStart, domainIncr, domainVal )
%
x=1+(domainVal-domainStart)./domainIncr;
idx=round(x);



