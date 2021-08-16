function  segdatDemo ()
% segdatDemo
% Demonstrate segdat and its functions
clear all

%% Generate a simple segdat oeject by concatonating two nldat objects with no overlap 
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001,'chanNames',{'X'});
XS=segdat(X);
Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',2, 'chanNames',{'Y'});
YS=segdat(Y);
S=cat(1, XS,YS);
plotHelper (X,Y,S);
%% reverse order of segments

Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
YS=segdat(Y);
S=cat(1,XS,YS);
plotHelper (X,Y,S);%% domain
%%
S=segdat;
X=nldat(ones(2000,1),'domainIncr',.001,'chanNames',{'X'},'domainStart',0);
XS=segdat(X);
Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
YS=segdat(Y);
S=cat(1, XS,YS);
plotHelper (X,Y,S);

%%  Basic operations
S2=S+1;
S3=S-2;
S=2*S;



%% Retrieve segments from segdat as nldat objects
disp(S)
disp (['Number of segments =' num2str(segCount(S))]); 
clf
X1=segGet(S,1); subplot (2,1,1); plot(X); title('Segment 1: X'); 
X2=segGet(S,2);subplot (2,1,2); plot(Y); title ('Segment 2: Y'); 
disp('Segment 1');
disp(X1);
disp('Segment 2');
disp(X2);

%% nldat returns a nldat object of all data with segments separate by nans

clf; 
disp('S');
disp(S)
disp('nldat(S)');
SN=nldat(S);
disp(SN)
plot(SN);title ('nldat(S)');
set(gca,'ylim',[0 2.1]);

%% convert nldat with nans into segdat

S2=segdat(SN);X
disp(S2);
plot(S2)
set(gca,'xlim',[0 2],'ylim',[0 2.1]);

%% Generate a segdat object by concatonating two nldat objects with overlap 
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001,'chanNames',{'X'});
XS=segdat(X);
Y=nldat(ones(1500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
xLim=[0 3];
YS=segdat(Y);
S=segdat.cat(1, XS,YS);
clf;
subplot (3,1,1);
plot (X);
set(gca,'xlim',xLim,'ylim',[0 2.1]);
title('X');
subplot (3,1,2);
plot(Y);
set(gca,'xlim',xLim,'ylim',[0 2.1]);
title('Y'); 
subplot (3,1,3);
plot(S); 
title ('segdat of X and Y ');
set(gca,'xlim',xLim,'ylim',[0 2.1]);
%% Polynomial operations
X=nldat(rand(5000,1),'domainIncr',.001,'chanNames',{'X'});
Y=X.^2;
Z=cat(2,X,Y);
S=segdat.randSeg(Z, 100, 200, 50);

P=polynom(Z,'polyType','power');
[R,V,ZP]=nlid_resid(P,S);

%% SSM models
x=rand(5000,1);
X=nldat(X,'domainIncr',.01);
Z=nlid_sim('L2',X);
ZS=segdat.randSeg(Z,200, 400, 50);
SS=ssm(ZS);
[R,V,ZP]=nlid_resid(SS,ZS);


%% NL BL Models 
x=randn(5000,1);
X=nldat(x,'domainIncr',.01);
Z=nlid_sim('N2L',X);
ZS=segdat.randSeg(Z,200,500,50);

NL=nlbl(ZS,'idMethod','subspace');
[R,V,P]=nlid_resid(NL,ZS);



%% Decimate
xd=decimate(S,2);clf; plot(xd);

%% Mulitple channel support 
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001,'chanNames',{'X'});
XS=segdat(X);
Y=nldat(ones(1500,1)*2,'domainIncr',.001,'domainStart',2., 'chanNames',{'Y'});
YS=segdat(Y);
S=segdat.cat(1, XS,YS);






end
function plotHelper (X,Y,S);
xLim=[0 3];
clf;
subplot (3,1,1);
plot (X);
set(gca,'xlim',xLim,'ylim',[0 2.1]);
title('X');
subplot (3,1,2);
plot(Y);
set(gca,'xlim',xLim,'ylim',[0 2.1]);
title('Y'); 
subplot (3,1,3);
plot(S); 
title ('segdat of X and Y ');
set(gca,'xlim',xLim,'ylim',[0 2.1]);
end