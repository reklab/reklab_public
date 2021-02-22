function  lnlDemo
%lnlDemo - demonstrate lnlbl

X=nldat(rand(5000,1),'domainIncr',.01);
X=X-mean(X);X=X*10; 
[Z,M,N]=nlid_sim('LNL',X);
figure(1); clf
plot(M);
title('Simulated LNL model');
%%%Iniialize paramters of model elements
IM=lnlbl;
set(IM,'initMethod','kernels'); 
i1=irf;set(i1,'nLags',51,'nSides',1);
p=polynom;set(p,'polyType','tcheb'); 
i2=irf; set(i2,'nLags',51,'nSides',1);
IM.elements= {i1 p i2};
IM=nlident(IM,Z);
figure(2);clf;
plot(IM); streamer ('Estimated LNL Model');


%% Predicted and observed outpus 
figure(3); clf

[R, V, yp]=nlid_resid(IM,Z,'plotflag',true)
disp(['%VAF=' num2str(double(V))]);
%% Estimate model with two sided IRfs 
IM=lnlbl;
i1=irf; set(i1,'nSides',2,'nLags',51)

p=polynom;
i2=irf; set(i2,'nSides',2,'nLags',51); 
IM.elements= {i1 p i2};
IM=nlident(IM,Z);


figure(2);clf;
plot(IM); streamer ('Estimated LNL Model');
figure(3);
[R, V, yp]=nlid_resid(IM,Z,'plotflag',true);
disp(V)
end

