function nlmtst(P)
% nltest function for pdf operater
disp('nlmtst for pdf');
figure(1)
x=randn(10000,1);
X=nldat(x);
P=pdf(X);
subplot (2,2,1);
plot(pdf(P,x,'pdfType','Density'));
subplot (2,2,2);
plot(pdf (P,x,'pdfType','Frequency'));
subplot (2,2,3);
plot(pdf (P,x,'pdfType','Probability'));
subplot(2,2,4)
plot(pdf (P,x,'pdfType','CumulativeProbability'));

figure(2);
subplot (2,1,1);
P=(pdf(X,'pdfType','Density','nBins',10,'binMin',-5,'binMax',5));
plot(P)



% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 
