function [ Y ]   = pcModel (x); 
% [ Y, M]   = pcModel (x); 
% Simulate a parallel cascade model consisting of a 
% LN system in parallel with a Weiner system. 
% First linear element
% x - input sampled at 1 kHZ
% Y - simulated output
% M - model 
% 
domainIncr=.001;
t=(0:domainIncr:.1)';
T=nldat(t,'domainIncr',domainIncr);
L1w=waveform('waveformType','irf2');
set(L1w,'gain',1,'damping',.8,'freq',450);
L1=nlsim(L1w,T);
d=cat(1, [0 0 0 0]',L1.dataSet);
L1=set(L1,'dataSet',d,'comment','First Linear Element');
%
x1=max(x,-.5);
x1=x1+ x1.^2;

y1=nlsim(L1,x1); 

L2w=waveform('waveformType','irf2');
set(L2w,'gain',-2,'damping',.4,'freq',100);
L2=nlsim(L2w,T);

y2=nlsim(L2,x);
y3= y2 + y2.^3;

%
Y=2*y1+y3;

