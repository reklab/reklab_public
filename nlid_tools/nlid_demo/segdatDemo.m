function  segdatDemo ()
% segdatDemo
% Demonstrate segdat and its functions
clear all
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001);
XS=segdat(X);
Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',1.5);
YS=segdat(Y);
S2=segCat(XS,YS);
plot(S2)

X1=setget(XS,1);

size(S)


plot(S