function nonlnearID_demo (model_type, x, noise_level);
% noninearID_demo - Demonstrate Object Oriented NLID ientification
% xin - input signal
% model_type - type of model to simulate; see nlid_sim for valid options;
%					 default value is 'LN2';
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt

% noise level - normalized noise power
nlevel=0;
ans='y';
pauseFlag=true;

if nargin==0 | isempty(option)
    opts = { 'LP1' 'HP1'  'L1' 'L2'  'LNRect' 'LN2' 'LN3' 'N3L' 'N2L' ...
        'LNL' 'PC' 'POLY' 'Static_Linear' 'Cuber' 'N3HP'};
    o=menu('Option',opts);
    modelType=opts{o};
end

%% default input signal is unit variance white noise
if nargin <2,
    inputAmp=1;
    inputAmp=input_d('Input Amplitude',inputAmp,1,100);
    x=randn (1000,1)*inputAmp;
    x=nldat(x,'domainIncr',.01);
end
% default noise level is 10%
if nargin < 3,
    noise_level=.1; %
    noise_level=input_d('noise level',noise_level,0,10);
end
pltflg=false;
% simulate noise-free response

[z,M]=nlid_sim (modelType,x,'plotFlg',pltflg);
%%
% generate noise corrupted output
%
y=double(z(:,2));
noise=randn(length(y),1);
scale=sqrt(noise_level*cov(y)/cov(noise));
yn= y + noise*scale;
z(:,2)=nldat(yn);
set(z,'chanNames',{'xIn' 'yOut'});
%%
%plot model
%pauseFlag=false;

figure (1); clf;
set(M,'comment','LNL Model'); 
plot (M); 
figure(2); plot(z,'nh',1);
streamer ('Input and output');
demo_pause

%%
% Linear impulse response function
%
nLags=32;
comment='Linear IRF Model Predictions';
disp(comment);
i=irf;
i=nlident(i,z,'nLags',nLags);
demo_plot(i,z, comment);

%% Fit a ln model
comment='LN Model';
disp(comment);
ln=lnbl;
ln{1,1}.nLags=nLags;
ln=nlident(ln,z);
demo_plot(ln,z,comment);



%%  Weiner series from ln estimate
W=wseries(ln);
set(W,'comment','Weiner series of LN model');
plot(W);
e=W.elements;
figure(2);
clf;
plot(e{3});
streamer ('Second order Weiner kernel from ln model');
demo_pause;


%%
% Fit an nl model
%
comment='NL Model';
disp(comment);
nl=nlbl;
nl{1,2}.nLags=nLags;
nl=nlbl(z)
disp(comment);
demo_plot(nl,z,'NL Model Predictions');

%% Weiner series from the NL model
W=wseries(nl);
set(W,'comment','Weiner series of NL model');
plot(W);

figure(2);clf
el=W.elements;
plot (el{3});
streamer ('Second order Weiner kernel from NL model');

demo_pause

%% Fit a nln model
comment='NLN Model';
disp(comment);
lnl=lnlbl;
nln{1,1}.nLags=nLags;
lnl=nlident(lnl,z);
set(lnl,'comment','Estimate LNL model'); 
demo_plot(lnl,z,comment)

%% Weiner series from the LNL model
W=wseries(lnl);
set(W,'comment','Weiner series of LNL model');
plot(W);

figure(2);clf
el=W.elements;
plot (el{3});
streamer ('Second order Weiner kernel from LNL model');

demo_pause

%% Estimate a parallel cascade model

comment=('Parallel Cascade model');
disp(comment);
pc = pcascade;
set(pc,'nPaths',5,'polyOrderMax',5,'nLags',nLags);
pc=nlident(pc,z);
set (pc,'comment','Parallel cascade model');
demo_plot(pc,z,comment);
%plot(wseries(pc));

%% Estimate the Weiner kernels using  Lee-Schetzen
comment='Wiener Kernel (LS) Model';
wk=wseries(z,'nLags',32,'method','LS','comment',comment);

figure(3);clf
el=wk.elements;
plot (el{3});
streamer ('Second order Weiner kernel from LS');
demo_plot(wk,z,comment');
%% Estimate Weiner kernels using the Toeplitz method
comment='Wiener Kernel (Toeplitz) Model';
disp(comment);
wt=wseries(z,'nLags',32,'method','Toeplitz','comment',comment);
figure(3);clf
el=wt.elements;
plot (el{3});
streamer ('Second order Weiner kernel from Toeplitz');
demo_plot(wt,z,' ')




    function demo_plot (model, z, comment)
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
        comment=get(model,'comment');
        commnet=comment(:)';
        
        set(zm,'comment',[comment ';Predictions  Vaf =' num2str(double(vf))]);
        plot(zm,'nv',3);
        %streamer(zm.comment);
        demo_pause
        return
    end



    function demo_pause
        if ~pauseFlag
            return
        end
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











