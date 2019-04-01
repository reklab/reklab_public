function polynomDemo 
% polydemo demonstrate for polynom objects
disp('polynomDemo - demnstrate and test use of polynom');

%% Test creation of polynominl fdisp('nlmtst for polynom objects');
pTypes = { 'tcheb' 'power' 'hermite'  'Bspline' 'laguerre'};

%% Dispay basis functions for different polynomials
nType=length(pTypes);
for iType=1:nType,
    curType=pTypes{iType};
    p=polynom('polyType',curType);
    B=basisfunction(p);
    figure(iType);
   
    titleStr=B.comment;
    set(B,'comment',''); 
     plot(B);
    streamer(titleStr,.90); 
end

%% Fit powet,tcheb and hermite to a cubic polynomial 

nx=3; ny=2;
F=1;
% Generate polynomial data set
%
x= (-10:.01:10)';
y= 10 + 5*x +2*x.^2 + .1*x.^3;
noise=randn(length(x),1)*10;
yNoise=y+noise;
z=cat(2,x,yNoise);
z=nldat(z,'domainIncr',.01);
%z=nlid_sim('poly');
subplot (nx,ny,1);
plot (z,'plotmode','xy');
title('data');
F=2;
for t=1:3,
    p=[];
    % Estimate polynomial
    figure(t);clf
    curType=pTypes{t};
    p=polynom(z,'polyType',curType,'polyOrder',10,'nInputs',1);
    subplot (nx,ny,1);
    plot(p);
    title ([ ' Estimated ' curType ])
    % show prediction
    subplot (nx,ny,2);
    
    yPre=nlsim(p,z(:,1));
    plot(z(:,2));
    h=line(yPre); set(h,'color','r');
    title (['Prediction for estimated: ' curType])
    % Create polynomial apriori
    ;
    % test a priori generation and simulation
    p1=polynom('polyType',curType,'polyCoef',p.polyCoef);
    set(p1,'polyMean',mean(x),'polyStd',std(x),'polyRange',[ min(x) max(x)]');
    subplot (nx,ny,3); plot(p1);
        title ([ ' a priori ' curType ])
    % prediction of a apriori polynomal
    subplot (nx,ny,4);
    yPre1=nlsim(p1,z(:,1));
    plot(z(:,2));
    h=line(yPre1); set(h,'color','r');
     title (['Prediction for apriori: ' curType])
    % Derivative
    subplot(nx,ny,5);
    plot (ddx(p));
    title('derivative'); 
    streamer(curType,.9); 
end

%%

for t=4:5,
    p=[];
    % Estimate polynomial
    figure(t);clf
    curType=pTypes{t};
    p=polynom(z,'polyType',curType,'polyOrder',10,'nInputs',1);
    subplot (nx,ny,1);
    plot(p);
    title ([ ' Estimated ' curType ])
    % show prediction
    subplot (nx,ny,2);
    
    yPre=nlsim(p,z(:,1));
    plot(z(:,2));
    h=line(yPre); set(h,'color','r');
    title (['Prediction for estimated: ' curType])
    % Create polynomial apriori
    ;
    % test a priori generation and simulation
    p1=polynom('polyType',curType,'polyCoef',p.polyCoef);
    set(p1,'polyMean',mean(x),'polyStd',std(x),'polyRange',[ min(x) max(x)]');
    subplot (nx,ny,3); plot(p1);
        title ([ ' a priori ' curType ])
    % prediction of a apriori polynomal
    subplot (nx,ny,4);
    yPre1=nlsim(p1,z(:,1));
    plot(z(:,2));
    h=line(yPre1); set(h,'color','r');
     title (['Prediction for apriori: ' curType])
    % Derivative
    subplot(nx,ny,5);
    plot (ddx(p));
    title('derivative'); 
    streamer(curType,.9); 
end
%%
return


% Test of double
pd = double(p);

% Hessian

% h= hessian(p,z);

% Jacobian
% j = jacobian(p,z);
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

%% B_splines

disp ('B-splines');
gain=zeros(200,1);
gain=cat(1,gain,linspace(0,5,600)');
gain=cat(1,gain, zeros(200,1));
G=nldat(gain,'domainIncr',.001);
BF=polynom;
set(BF,'polyType','Bspline','polyOrder',10);
BFnew=nlident(BF,G)
D=G;
set(D,'dataSet',domain(G))
GP=nlsim(BFnew,D);
figure(4);
clf
plot(G)
h=line(GP);set(h,'color','r');
title('B-spline demo');
legend ('Signal','Spline fit'); 




% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt
