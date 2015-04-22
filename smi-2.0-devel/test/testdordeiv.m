% example taken from smidemo5
% The plant is a mechanical system consisting of 2 rotating discs
% and an integrator. Its transfer function (order = 5) is 
Pnum = [0 0.98 12.99 18.59 3.30 -0.02]*1e-3; 
Pden = [1 -4.4 8.09 -7.83 4 -0.86]; 

% A controller has been designed based on H_infinity method
Cnum = [0.61 -2.03 2.76 -1.83 0.49];
Cden = [1 -2.65 3.11 -1.75 0.39];

% The process noise is shaped by the following filter 
Fnum = [0 2.89 11.13 2.74]*0.01;
Fden = [1 -2.7 2.61 -0.9];

plant=tf(Pnum,Pden,-1);
controler=tf(Cnum,Cden,-1);
noisefilter=tf(Fnum,Fden,-1);
interconnection=ss([],[],[],[1 0 1; 0 1 0],-1); 
combinedsystem1=parallel(plant,noisefilter,[],[],1,1);
combinedsystem=parallel(combinedsystem1,interconnection,1,1,1,2);
closedloopsystem=feedback(combinedsystem,controler,2,1,-1);

N = 1200; 
r = randn(N,1);      % external input 
v= 0.01*randn(N,1);  %output measurement noise
f= 0.01*randn(N,1);  %input measurement noise
w=(1/3)*randn(N,1);  % proces noise
Ucl=[w,r,v,f];
Ycl = lsim(closedloopsystem,Ucl);
u = Ycl(:,2);
y = Ycl(:,1); 

[S,R]=dordeiv(u,y,[],20);
[S,R]=dordeiv(u,y,r,20);
[S,R2]=dordeiv(u,y,r,20,R);

S



