function  lnlDemo
%lnlDemo - demonstrate lnlbl

X=nldat(rand(5000,1),'domainIncr',.01);
X=X-mean(X);
[Z,M,N]=nlid_sim('LNL',X);
figure(1); clf
plot(M);
title('Simulated LNL model');
% Estimate model
IM=lnlbl(Z);
figure(2);clf;
plot(IM); streamer ('Estimated LNL Model');
%% Predicted and ovbserved outpus 
figure(3); clf

[R, V, yp]=nlid_resid(IM,Z,'plotflag',true)


end

