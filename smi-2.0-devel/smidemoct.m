 n=4; m=1; l=2;
As= [-0.49  0.00005  -4.8   0
    0    -0.015   -14   -32.2 
    1    -0.00019 -0.84   0
    1     0        0      0];
Bs= [-8.7;-1.1;-0.11;0];
Fs= [-4.8;-14;-0.84 ;0];
Cs= [1 0 0 0;
    0 1 0 0];
Ds=zeros(2,1); 
Te=500;
Ts=0.05;
t=(0:Ts:Te-Ts)';
N=size(t,1);
omega=[0.5;1;2;5;10];
phase =[0.0245;-1.3488;1.7848;-0.9223;0.0061];
u=zeros(N,1);y=zeros(N,2);
for i=1:5
  [ui,yi]=aircraft(omega(i),t,phase(i));
  u=u+ui;y=y+yi;
end


u2=zeros(N,1);
for i=1:5
  u2=u2+sin(omega(i)*t+phase(i));
end
np=randn(N,1)*0.2;% alpha_g
nm=randn(N,2)*diag([0.2,0.05]);
y2=lsim(As,[Bs Fs],Cs,[Ds zeros(2,1),[u np],t)+nm;
  

yd=lsim(As,Fs,Cs,Ds,np,t)+nm;

[Sc,Rc]=cordpo(u2,y2,Ts,4,30);
[Ac,Cc]=cmodpo(R1,2);
[Bc,Dc]=cac2bd(Ac,Cc,u,y,Ts);


