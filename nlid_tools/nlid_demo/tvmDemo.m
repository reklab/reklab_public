% tvmDemo - demonstrate tvm propoerties
%% TV IRF
%Generate an TV IRF
gain=ones(200,1);
gain=cat(1,gain,linspace(1,5,600)');
gain=cat(1,gain, ones(200,1));
plot(gain)
tvI={};
 for i=1:1000,
     tvTemp=irf2(irf,'g',gain(i));
     tvI{i}=tvTemp;
 end
 TVM=tvm;
 
set(TVM,'elements',tvI,'tvStart',0,'tvIncr',.01);
figure(1);
plot(TVM);

x=randn(1000,1,200);
X=nldat(x,'domainIncr',.01);

Y=nlsim(TVM,X);
Z=cat(2,X,Y); 

I=irf;
set(I,'nLags',50);
tvIdent=nlident(tvm,Z,I);

simY=nlsim(tvIdent,X);


; 
   