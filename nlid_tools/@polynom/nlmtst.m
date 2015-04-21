function [p,z] = nltst(p);
% polynom/nltst for polynom objects

%% Test creation of polynominl fdisp('nlmtst for polynom objects');
ptypes = { 'tcheb' 'power' 'hermite' };
nx=3; ny=2;
F=1;
% Generate polnomial data set
%
x= (-10:.01:10)';
y= 10 + 5*x +2*x.^2 + .1*x.^3;
noise=randn(length(x),1)*10;
yNoise=y+noise;
z=cat(2,x,yNoise);
z=nldat(z);
%z=nlid_sim('poly');
subplot (nx,ny,1);
plot (z,'plotmode','xy');
title('data');
F=2;
for t=1:length(ptypes),
    p=[];
    % Estimate polynomial
    figure(t);clf
    p=polynom(z,'polyType',ptypes{t},'polyOrder',3,'nInputs',1);
    subplot (nx,ny,1);
    plot(p);
    title ([ ' Estimated ' ptypes{t} ])
    % show prediction
    subplot (nx,ny,2);
    
    yPre=nlsim(p,z(:,1));
    plot(z(:,2));
    h=line(yPre); set(h,'color','r');
    title (['Prediction for estimated: ' ptypes(t)])
    % Create polynomial apriori
    ;
    % test a priori generation and simulaation
    p1=polynom('polyType',ptypes{t},'polyCoef',p.polyCoef);
    set(p1,'polyMean',mean(x),'polyStd',std(x),'polyRange',[ min(x) max(x)]');
    subplot (nx,ny,3); plot(p1);
        title ([ ' a priori ' ptypes{t} ])
    % prediction of a apriori polynomal
    subplot (nx,ny,4);
    yPre1=nlsim(p1,z(:,1));
    plot(z(:,2));
    h=line(yPre1); set(h,'color','r');
     title (['Prediction for apriori: ' ptypes(t)])
    % Derivative
    subplot(nx,ny,5);
    plot (ddx(p));
    title('derivative'); 
end

%%

% Test of double
pd = double(p);

% Hessian

h= hessian(p,Z);

% Jacobian
j = jacobian(p,Z);
return
% nlid_resid(p,z);

%
% Now try two input third order
x1=x(20:-1:1);
y2 = 10 + 1.1*x + 1.2*x1 + 2.1*x.^2 + 2.2*x.*x1 + 2.3*x1.^2 + ...
    + 3.2*x.^2.*x1 + 3.3*x.*x1.^2
z=cat(2,x',x1',y2');
p=polynom(z,'order',3,'NInputs',2,'type',T);
figure(5);
F=F+1;
clf;
plot (p);
figure(6);
clf;
nlid_resid(p,z);


% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt
