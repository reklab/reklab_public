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
S=cat(XS,YS);
plotHelper (X,Y,S);
%% reverse order segmens

Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
YS=segdat(Y);
S=cat(XS,YS);
plotHelper (X,Y,S);%% domain
%%
S=segdat;
X=nldat(ones(2000,1),'domainIncr',.001,'chanNames',{'X'},'domainStart',0);
XS=segdat(X);
Y=nldat(ones(500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
YS=segdat(Y);
S=cat(XS,YS);
plotHelper (X,Y,S);

%% 


%% Retrieve segments from segdat as nldat objects
disp(S)
disp (['Number of segments =' num2str(segCount(S2))]); 
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
plot(SN);title ('nldat(S2)');
set(gca,'xlim',xLim,'ylim',[0 2.1]);

%% convert nldat with nans into segdat

S2=segdat(SN);X
disp(S2);
set(gca,'xlim',[0 2],'ylim',[0 2.1]);

%% Generate a segdat oeject by concatonating two nldat objects with overlap 
S=segdat;
X=nldat(ones(1000,1),'domainIncr',.001,'chanNames',{'X'});
XS=segdat(X);
Y=nldat(ones(1500,1)*2,'domainIncr',.001,'domainStart',.5, 'chanNames',{'Y'});
YS=segdat(Y);
S=segCat(XS,YS);
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

%% Mulitple channel support 
NOT YET IMPLEMENTED





%% Decimate
xd=decimate(S2,2);clf; plot(xd);

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