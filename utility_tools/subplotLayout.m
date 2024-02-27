function [nH,nV]=subplotLayout (nSig)
% [nH,nV]=subplotLayout (nSig)
% Select subplot layout for mulutiple plots
if nSig<=4
    nH=1;
    nV=nSig;
elseif nSig<=6
    nH=2;
    nV=3;
elseif nSig<=8
    nH=2;
    nV=4;
elseif nSig<=12
    nH=3;
    nV=4;
elseif nSig<=16
    nH=4;
    nV=4;
else
    nH=ceil(sqrt(nSig));
    nV=ceil(nSig/nH);
end
end