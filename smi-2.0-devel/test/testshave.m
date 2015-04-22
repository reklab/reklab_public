y=rand(1000,1);
n1=rand(10,1);
n2=rand(10,1);

y(100:109)=10*n1;
y(500:509)=5*n2;
ys=shave(y)
ys=shave(y,2,0.05)
subplot(2,1,1)
plot(y)
title('peaky signal')
subplot(2,1,2)
plot(ys)
title('shaved signal')

