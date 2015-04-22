function x=linAxis (domainStart, domainIncr, domainEnd)
% Generrate a nlda tobject correspindiong to linear axis
x=domainStart:domainIncr:domainEnd;
x=nldat(x','domainIncr',domainIncr,'domainStart',domainStart,'chanNames',{'linAxis'});
return

