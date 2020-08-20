function  segdatDemo ()
% segdatDemo
% Demonstrate segdat and its functions
clear all
%% Generate a simple segdat obeject by concatonating two nldat objects
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001);
XS=segdat(X);
Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',1.5);
YS=segdat(Y);
S2=segCat(XS,YS);
clf;
plot(S2); set(gca,'ylim',[0 3]); title ('segdat of X and Y ');

%% Retrieve segments as nldat objects

clf
X1=segGet(S2,1); subplot (2,1,1); plot(X); title('X'); 
X2=segGet(S2,2);subplot (2,1,2); plot(Y); title ('Y'); 

%% nldat returns a nldat object of all data. 

clf; plot(nldat(S2));

%% Mulitple channel support 
X2=cat(2,X,-X);
Y2=cat(2,Y,-Y);
Z2=segdat(X2);
Z2=cat(1,Z2,Y2);
segPlot(Z2);

Z3=cat(2,XS,-XS);



%% domain
d=domain(Z2);
clf; plot (d)

size(S)


xd=decimate(S2,2);clf; plot(xd);

