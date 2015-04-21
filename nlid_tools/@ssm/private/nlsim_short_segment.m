function outp = nlsim_short_segment (ps,z)
%This function estimates initial condition of a linear system and predicts
%model output using the estimated initial conditions
%
%%
assign(ps);
ts = get(z,'domainIncr');
if ~(ts == domainIncr)
    error('Input-output and system sampling rate mismatch')
end
in_onsetPointer = get(z,'onsetPointer');
onsetPointer = in_onsetPointer (:,2);
in_onsetPointer = in_onsetPointer (:,1);
in_segLength = get(z,'segLength');
segLength = in_segLength (:,2);
in_segLength = in_segLength (:,1);
if ~( isequal(onsetPointer,in_onsetPointer) &&  isequal(in_segLength,segLength))
    error('The input and output onset pointer and length must be equal..')
end
data = get(z,'dataSet');
input = data(:,1);
output = data(:,2);
endpointer = onsetPointer + segLength - 1;
%extracting input-output data from segdat
in = zeros(sum(segLength),1);
out = zeros(sum(segLength),1);
pointer = 1;
switch_time = zeros(length(endpointer)-1,1);
for i = 1 : length(endpointer)
    in_temp = input(onsetPointer(i):endpointer(i));
    in_temp = delay(nldat(in_temp,'domainIncr',ts),nDelayInput);
    in_temp = get(in_temp,'dataSet');
    in(pointer:pointer+segLength(i)-1) = in_temp;
    out(pointer:pointer+segLength(i)-1) =output(onsetPointer(i):endpointer(i));
    pointer = pointer + segLength(i);
    switch_time(i) = pointer;
end
out = out - mean(out);
in = in - mean(in);
switch_time = [1;switch_time];
nsamp = length(in);
p = length(segLength);
m = size(A,1);
%Defining regressor matrices
%Gamma is the regressor for the initial conditions
Gamma_total = zeros(size(in,1),p*m);
Phi_total = zeros(size(in,1),(m+1));
max_interval = max(segLength);
Gamma_nominal = zeros(max_interval,m);
Gamma_nominal(1,:) = C;
An =A;
for i = 1:floor(log(nsamp)/log(2))
    Gamma_nominal(2^(i-1)+1:2^i,:) = Gamma_nominal(1:2^(i-1),:)*An;
    An = An * An;
end
Gamma_nominal(2^i+1:nsamp,:) = Gamma_nominal(1:nsamp-2^i,:) * An;
for i = 1 : p
    Gamma_total(switch_time(i):switch_time(i+1)-1,(i-1)*m+1:i*m) = Gamma_nominal(1:segLength(i),:);
    Phi = BD_omega_regressor(in(switch_time(i):switch_time(i+1)-1,:),A,C);
    Phi_total(switch_time(i):switch_time(i+1)-1,:) = Phi;
end
initial = lscov(Gamma_total,out-Phi_total*[B;D]);
initial = reshape(initial,m,length(initial)/m);
output_predicted = zeros(size(out));
for i = 1 : p
    output_predicted(switch_time(i):switch_time(i+1)-1) = dlsim(A,B,C,D,in(switch_time(i):switch_time(i+1)-1,:),initial(:,i));
end
output_predicted = output_predicted - mean(output_predicted);
outp = segdat;
chanName = get(z,'chanNames');
chanUnit = get(z,'chanUnits');
set(outp,'dataSet',output_predicted,'chanNames','Predicted output','chanUnits',chanUnit,'chanNames',chanName{2},'domainIncr',ts,'onsetPointer',switch_time(1:end-1),'segLength',segLength);
% figure
% plot(output_predicted,'r');
% hold on
% plot(out)
% output_predicted = nldat(output_predicted,'domainIncr',ts);
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