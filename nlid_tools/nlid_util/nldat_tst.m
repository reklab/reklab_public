function  nldat_tst( X )
%UNTITLED Test opweration of nldat objects 
%   Detailed explanation goes here
x=randn(1000,1);
X=nldat(x, 'comment','Test comment');
set (X,'domainIncr',.001,'chanNames',{'One'}, 'chanUnits',{'Volts'});
set(X,'domainName','Time','domainStart',0); 
plot (X);
plot (chop(X,.1))
plot (cumsum(X))
plot (ddt(X))
plot(decimate(X,10))
plot (detrend(X));
disp(X);
plot (domain(X))
plot (double(X))
% emean
plot (ext(X,.1,.5));
plot (extract(X,100,500));
f=fft(X);
plot(abs(f));
plot (angle(f));
s= get(X);

plot(hwrect(X));
plot(isnan(X));
% iddata
length(X);
max(X);
min(X);

mean(X);
std(X);

plot (X+X);
plot (X+5);
plot (X-X); 
plot(X-5)
plot (10*X)
plot (X./10);
Y=reshape(X,50,2,10);
plot(Y)
plot(emean(Y));
stairs(X);
plot (X(100:200,1,1))
plot (Y(:,1,1))
X1=X;
X1(100:200)=ones(101,1); plot(X1)

var(X);
Z=cat(2,X,X);
vaf(X,X);
end