function nlid_demo (model_type, x, noise_level);
% NLID_DEMO - Demonstrate Object Oriented NLID ientification
% xin - input signal
% model_type - type of model to simulate; see nlid_sim for valid options;
%					 default value is 'LN2';
%


% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt
% Seocnd trst comment
% noise level - normlaized noise power
nlevel=0;
ans='y';
if nargin == 0,
    model_type= 'LN2';
end

%% default input signal is unit variance white noise
if nargin <2,
    x=randn (1000,1);
    x=nld
    at(x,'domainIncr',.01);
end
% default noise level is 10%
if nargin < 3,
    noise_level=.1; %
end
pltflg=false;
% simulate noise-free response

[z,M]=nlid_sim (model_type,x,'plotFlg',pltflg);
%
% generate noisey output measure
%
y=double(z(:,2));
noise=randn(length(y),1);
scale=sqrt(nlevel*cov(y)/cov(noise));
yn= y + noise*scale;
z(:,2)=nldat(yn);
%
% plot model
figure (1); clf;
plot (M); title([ model_type ' Model']);
demo_pause
figure(1);  plot(z,'nh',1);
demo_pause

%
% Linear impulse response function
%
nLags=32;
comment='Linear IRF Model';
disp(comment);
i=irf;
i=nlident(i,z,'nLags',nLags);
demo_plot(i,z, comment);
%
% Fit a ln model
%
comment='LN Model';
disp(comment);
ln=lnbl;
ln=nlident(ln,z, 'nLags',nLags);
demo_plot(ln,z,comment);
W=wseries(ln);
plot(W);
e=W.elements;
figure(2);
clf;
plot(e{3});
demo_pause;


%
% nl model
%
comment='NL Model';
disp(comment);
nl=nlbl(z,'nLags',nLags);
disp(comment);
demo_plot(nl,z,comment);
W=wseries(nl);
disp('Wiener series of NL model');
plot(W);
figure(2);clf
el=W.elements;
plot (el{3});
demo_pause
%
%
% pcascade
%
comment=('Parallel Cascade model');
disp(comment);
pc = pcascade;
set(pc,'nPaths',5,'polyOrderMax',5,'nLags',nLags);
pc=nlident(pc,z);
set (pc,'comment','Parallel cascade model');
demo_plot(pc,z,comment);
%plot(wseries(pc));

% wkernel = Lee-Schetzen
comment='Wiener Kernel (LS) Model';
disp(comment);
wk=wseries(z,'nLags',32,'method','LS');
demo_plot(wk,z,comment);
%
% wkernel - Toeplitz
comment='Wiener Kernel (Toeplitz) Model';
disp(comment);
wt=wseries(z,'nLags',32,'method','Toeplitz');
demo_plot(wt,z,' ')
% Voltera kernel (Fast Orthogonal)
% Not working properly
% if input_l('Volterra kernel (slow!)','N');,
%     comment='Voltera kernel (Fast Orthogonal)';
%     disp(comment);
%     vk=vseries(z,'nLags',20);
%     set(vk,'comment',comment)
%     demo_plot(vk,z,' ');
% end
% %
% % Wiener kernel (WB );
% comment='Wiener Bose (Laguerre) Model';
% disp(comment);
% wb=wbose(z,'Nlags',32,'Nfilt',10,'OrderMax',2);
% demo_plot(wb,z,' ');
% %
%
%    disp('Done');
%

function demo_plot (model, z, comment)
 figure(2); 
u=z(:,1);
y=z(:,2)-mean(y);
yp=nlsim(model,u);
yp=extract(yp,100);
y=extract(y,100);
yr=y-yp;
vf=vaf(y,yp);
zm=cat(2,y,yp);
zm=cat(2,zm,yr);
set(zm,'chanNames',{'Observed' 'Predicted' 'Residuals'});
resid=y-yp;
figure(1);
clf;subplot (1,1,1);
plot(model);


figure(2);clf
set(zm,'comment',[comment '; Vaf =' num2str(double(vf))]);
plot(zm,'nv',3);
demo_pause
return
end



function demo_pause
figure(1);
p1=get(gcf,'position');
p2=p1;
p2(1)=p1(3)-100;
p2(2)=p1(4)-25;
p2(3)= 80;
p2(4)=20;
h = uicontrol('Position',p2,'String','Continue',...
    'Callback','uiresume(gcbf)', 'backgroundcolor','red');
uiwait(gcf);clf
end
end











