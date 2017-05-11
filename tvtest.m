I={};

for i=1:200,
    if i<100,
        G=2;
    else
        G = 2 - (100-i)/100;
    end
    I{i}=irf2(irf,'g',G);
end
ITV=tvm;
set (ITV,'elements',I,'tvStart',1.01,'tvIncr',.01);
x=randn(1000,100);
x=reshape(x,1000,1,100);

X=nldat(x,'domainIncr',.01);

 ySimTV = nlsim(ITV, X);
  Z=cat(2,X, ySimTV);
  iTVest=tvm (Z,irf( 'nSides',1,'nLags',101));
   figNum=figNum+1; figure(figNum);

  plot(nldat(iTVest(:,:,100:end),'realizationMode','mesh'))
  