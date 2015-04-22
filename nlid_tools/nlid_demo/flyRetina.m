function [ Z,M ] = flyRetina( X )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
domainIncr=.001;
t=(0:domainIncr:.025)';
T=nldat(t,'domainIncr',domainIncr);

% First linear element
L1w=waveform('waveformType','irf2');
set(L1w,'gain',1,'damping',.8,'freq',450);
L1=nlsim(L1w,T);
d=cat(1, [0 0 0 0]',L1.dataSet);
L1=set(L1,'dataSet',d,'comment','First Linear Element');

% Second linear element

L2w=waveform('waveformType','irf2');
set(L2w,'gain',10,'damping',.5,'freq',800);
L2=nlsim(L2w,T);
d=cat(1, [0 0 0]',L2.dataSet);
L2=set(L2,'dataSet',d,'comment','Second Linear Element')
;
% first nonlinearity
x=[ -5 -15; -2.5 -7.5; 0 0  ; 2.5 2 ; 5 4]
z=nldat(x);
P1=polynom(z,'polyType','power','polyOrderMax',2','comment','First static Nonlinearity');

% second nonlinearity
x=[ -30 20; -20 0; -10 -10 ; 0 -5 ; 5 0; 10 6.5 ; 20 10 ; 30 15; 40 18]
z=nldat(x);
P2=polynom(z,'polyType','power','polyOrderMax',4,'comment','Second static nonlineairty');

% Define model
el={ L1 P1 L2 P2};
M=nlm;
set (M,'elements',el)
% Simulate response

Y=nlsim(M,X);
Z=cat(2,X,Y);









end

