function output_predicted = nlsim_hammer_subspace_short_segment (sys,z)
%This function estimates output of a hammerstein system based on short segments of input
%and output
%system= nlsim_hammer_subspace_short_segment (z)
%This routine is based on the following work:
%Kian Jalaleddini, Ferryl Alley, Robert E Kearney, "Identification of
%Hammerstein Systems from Short Segments of Data: Application to Stretch
%Reflex Identification", In proceeding of Sysid 2012, 16th IFAC Symposium
%on System Identification, Brussels, 2012.
% Author: Kian Jalaleddini
%%
elements = get(sys,'elements');
nonlinear = elements{1};
S = elements{2};
AT = get(S,'A');
BT = get(S,'B');
CT = get(S,'C');
DT = get(S,'D');
ts = get(z,'domainIncr');
delay = get(S,'nDelayInput');
m = size(AT,1);
z = z(:);
input = z(1:length(z)/3);
output_noisy = z(length(z)/3+1:2*length(z)/3);
input = input - mean(input);
output_noisy = output_noisy - mean(output_noisy);
switch_time = z(2*length(z)/3+1:3*length(z)/3);
switch_time = get(switch_time,'dataSet');
%Signal pre-processing
switch_time = switch_time(:);
zero_index = find(switch_time==0);
switch_time = switch_time(1:zero_index(1)-1);
switch_time = [1;switch_time];
p = length(switch_time);
switch_time = [switch_time;length(input)+1];
output_noisy = get(output_noisy,'dataSet');
input = nlsim(nonlinear, input);
input  = get(input ,'dataSet');
input = del(input,delay/ts);
interval = switch_time(2:end) - switch_time(1:end-1);
nsamp = length(input);
%Defining regressor matrices
%Gamma is the regressor for the initial conditions
Gamma_total = zeros(size(input,1),p*m);
Phi_total = zeros(size(input,1),(m+1));
max_interval = max(interval);
Gamma_nominal = zeros(max_interval,m);
Gamma_nominal(1,:) = CT;
An =AT;
for i = 1:floor(log(nsamp)/log(2))
    Gamma_nominal(2^(i-1)+1:2^i,:) = Gamma_nominal(1:2^(i-1),:)*An;
    An = An * An;
end
Gamma_nominal(2^i+1:nsamp,:) = Gamma_nominal(1:nsamp-2^i,:) * An;
for i = 1 : p
    Gamma_total(switch_time(i):switch_time(i+1)-1,(i-1)*m+1:i*m) = Gamma_nominal(1:interval(i),:);
    Phi = BD_omega_regressor(input(switch_time(i):switch_time(i+1)-1,:),AT,CT);
    Phi_total(switch_time(i):switch_time(i+1)-1,:) = Phi;
end
initial = lscov(Gamma_total,output_noisy-Phi_total*[BT;DT]);
initial = reshape(initial,m,length(initial)/m);
output_predicted = zeros(size(output_noisy));
for i = 1 : p
    output_predicted(switch_time(i):switch_time(i+1)-1) = dlsim(AT,BT,CT,DT,input(switch_time(i):switch_time(i+1)-1,:),initial(:,i));
end
output_predicted = output_predicted - mean(output_predicted);
output_predicted = nldat(output_predicted,'domainIncr',ts);
end
function Phi = BD_omega_regressor(u,A,C)
b=zeros(size(u,1),size(A,1)*(size(u,2)));
e=eye(size(A,1));
for j=1:size(u,2)
    for i=1:size(A,1)
        x=ltitr(A,e(:,i),u(:,j));
        yij=C*x';
        b(:,(j-1)*size(A,1)+i)=yij(:);
    end
end
bnew = zeros(size(b,1),size(b,2)+size(u,2));
k1=1;
k2=1;
for i=1:size(bnew,2)
    if mod(i,size(A,1)+1)==0
        bnew(:,i) = u(:,k1);
        k1 = k1+1;
    else
        bnew(:,i) = b(:,k2);
        k2 = k2+1;
    end
end
Phi = bnew;
end