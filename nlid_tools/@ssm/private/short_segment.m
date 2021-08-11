function system_ss = short_segment (z,ps)
%This function estimates a linear system based on short segments of input
%and output
%system_ss= linear_short_segment (z,ps)
%z must be a segdat object
%system_ss is an ssm object
% Author: Kian Jalaleddini
% Date: February 14, 2014 Ver0.1s
% Date : February 24, 2014 Ver0.2
%%
assign(ps);
condition = 1;
% changes to new segdat format where all channels must have same properties
ts = get(z,'domainIncr');
in_onsetPointer = get(z,'onsetPointer');
onsetPointer = in_onsetPointer ;% (:,2);
in_onsetPointer = in_onsetPointer; % (:,1);
in_segLength = get(z,'segLength');
segLength = in_segLength; %(:,2);
in_segLength = in_segLength; % (:,1);
% if ~( isequal(onsetPointer,in_onsetPointer) &&  isequal(in_segLength,segLength))
%     error('The input and output onset pointer and length must be equal..')
% end
data = get(z,'dataSet');
input = data(:,1);
output = data(:,2);
%Ensure each segment has enough number of smaples o.w. remove that segment
N = segLength - 2 * hankleSize + 1;
if length(find(N<1))>1
    warning(['Removing ',num2str(length(find(N<1))),' very short segments'])
end

onsetPointer(N<1) = [];
segLength(N<1) = [];
N(N<1) = [];
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
switch_time = [1;switch_time(:)];
nsamp = length(in);
p = length(segLength);
%Ensure enough number of samples is available
if nsamp>2*hankleSize*p-p+5*hankleSize+1
%Identify AT and CT 
%Define extended Hankle matrices for future-past input and future output
    Yf_tot = zeros(sum(N),hankleSize);
    Uf_tot = zeros(sum(N),hankleSize);
    Up_tot = zeros(sum(N),hankleSize);
    for i = 1 : p
        Uf = zeros(N(i), hankleSize); 
        Up = zeros(N(i), hankleSize); 
        Yf = zeros(N(i), hankleSize);
        input_segment = in(switch_time(i):switch_time(i+1)-1,:);
        output_segment = out(switch_time(i):switch_time(i+1)-1,:);
        for k = (1:hankleSize)
          Up(:,(k-1) +1:k ) = input_segment(k:N(i)+k-1,:); 
          Uf(:,(k-1) +1:k ) = input_segment(hankleSize+k:N(i)+hankleSize+k-1,:); 
          Yf(:,(k-1) +1:k ) = output_segment(hankleSize+k:N(i)+hankleSize+k-1,:); 
        end
        Yf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Yf;
        Uf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Uf;
        Up_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Up;
    end
    data_matrix = [Uf_tot Up_tot Yf_tot];
    [~ , R] = qr(data_matrix);
    L = R';
    L32 = L(2*hankleSize+1:3*hankleSize,hankleSize+1:2*hankleSize);
    [Un,Sn,~] = svd(L32); 
    Sn = diag(Sn); 
    Sn = Sn(1:hankleSize); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hankleSize);
%select linear system order    
    if strcmp(orderSelect,'manual')==1
        m = orderselect(Sn,'manual');
    elseif strcmp(orderSelect,'largest-gap')==1
        m = orderselect(Sn,'largest-gap');
    else
        error('orderselection must be set to either on or off')
    end
    if m==0 
        system_ss = ssm;
        condition = 0;
    end
    [AT , CT] = destac(R,m);
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Identified system is unstable.')
        warning('Attempt to identify the system failed.')
        system_ss = ssm;
        condition = 0;
    end
else
    warning('Not enough samples available for identification')
    warning('Attempt to identify the system failed.')
    system_ss = ssm;
    condition = 0;
end
if condition>0
%Defining regressor matrices
%Gamma_tot is the regressor for the initial conditions
%Phi_tot is the regressor for B & D state-space matrices
    Gamma_total = zeros(size(out,1),p*m);
    Phi_total = zeros(size(out,1),(m+1));
    max_interval = max(segLength);
    Gamma_nominal = zeros(max_interval,m);
    Gamma_nominal(1,:) = CT;
    An =AT;
    for i = 1:floor(log(nsamp)/log(2))
        Gamma_nominal(2^(i-1)+1:2^i,:) = Gamma_nominal(1:2^(i-1),:)*An;
        An = An * An;
    end
    Gamma_nominal(2^i+1:nsamp,:) = Gamma_nominal(1:nsamp-2^i,:) * An;
    for i = 1 : p
        Gamma_total(switch_time(i):switch_time(i+1)-1,(i-1)*m+1:i*m) = Gamma_nominal(1:segLength(i),:);
        Phi = BD_omega_regressor(in(switch_time(i):switch_time(i+1)-1,:),AT,CT);
        Phi_total(switch_time(i):switch_time(i+1)-1,:) = Phi;
    end
    Phi = [Gamma_total Phi_total];
    params = lscov(Phi,out);
    BT = params(p * m + 1:end-1);
    DT = params(end);
    initial = reshape(params(1 : p * m),m,length(params(1 : p * m))/m);
    system_ss = ssm;
    set(system_ss,'A',AT,'B',BT,'C',CT,'D',DT,'domainIncr',ts,'nDelayInput',nDelayInput);
    if displayFlag == 1
        output_predicted = zeros(size(out));
        for i = 1 : p
            output_predicted(switch_time(i):switch_time(i+1)-1) = dlsim(AT,BT,CT,DT,in(switch_time(i):switch_time(i+1)-1,:),initial(:,i));
        end
        output_predicted = output_predicted - mean(output_predicted);
        output_predicted = nldat(output_predicted,'domainIncr',ts);
        h = figure;
        for i =1 : p
            figure(floor((i-1)/4)+h)
            subplot(4,1,mod(i-1,4)+1)
            measured_data = out(switch_time(i):switch_time(i+1)-1);
            predicted_data = output_predicted(switch_time(i):switch_time(i+1)-1);
            predicted_data = predicted_data.dataSet;
            predicted_data = nldat(predicted_data,'domainIncr',ts);
            measured_data = nldat(measured_data,'domainIncr',ts);
            set(measured_data,'chanNames','Measured output');
            set(predicted_data,'chanNames','Predicted output');
            z_plot = cat(2,measured_data,predicted_data);
            set(z_plot,'chanUnits',get(z,'chanUnits'),'comment','Measured vs predicted output');
            plot(z_plot,'plotmode','super');
        end
    end
end
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