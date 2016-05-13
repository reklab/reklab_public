function system_nlbl = hammer_subspace_short_segment (z,ps)
%This function estimates a hammerstein system based on short segments of input
%and output
%system= hammer_subspace_short_segment (z)
%This routine is based on the following work:
%Kian Jalaleddini, Ferryl Alley, Robert E Kearney, "Identification of
%Hammerstein Systems from Short Segments of Data: Application to Stretch
%Reflex Identification", In proceeding of Sysid 2012, 16th IFAC Symposium
%on System Identification, Brussels, 2012.
% options={{'hankle_size' 20 'Size of hankle matrix'} ...
%          {'order' 8 'maximum order for nonlinearity'} ...
%          {'delay' 0 'Delay added to the input'} ...
%          {'orderDetectionMethod','preset'}...
%   { orderLe - order of LE 
%          {'threshold' 10^(-5)}...
%          {'plot_mode','No'}...
%      };
% if arg_parse(options,varargin);
%      return
%  end
%%
plotMode = 0;
assign(ps);
condition = 1;
ts = get(z,'domainIncr');
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
switch_time = [1;switch_time];
nsamp = length(in);
p = length(segLength);
if nsamp>2*hankleSize*p-p+2*maxOrderNLE*hankleSize+3*hankleSize+1
%Construct the input signal
    avg = (max(in) + min(in)) / 2;
    rng = max(in) - min(in);
    un = (in - avg) * 2 / rng;
    u = multi_tcheb(un,maxOrderNLE - 1);
%Identify AT and CT
    Yf_tot = zeros(sum(N),hankleSize);
    Uf_tot = zeros(sum(N),maxOrderNLE*hankleSize);
    Up_tot = zeros(sum(N),maxOrderNLE*hankleSize);
    for i = 1 : p
        Uf = zeros(N(i), maxOrderNLE*hankleSize); 
        Up = zeros(N(i), maxOrderNLE*hankleSize); 
        Yf = zeros(N(i), hankleSize);
        input_segment = u(switch_time(i):switch_time(i+1)-1,:);
        output_segment = out(switch_time(i):switch_time(i+1)-1,:);
        for k = (1:hankleSize)
          Up(:,(k-1) * (maxOrderNLE)+1:k * (maxOrderNLE)) = input_segment(k:N(i)+k-1,:); 
          Uf(:,(k-1) * (maxOrderNLE)+1:k * (maxOrderNLE)) = input_segment(hankleSize+k:N(i)+hankleSize+k-1,:); 
          Yf(:,(k-1) * 1+1:k * 1) = output_segment(hankleSize+k:N(i)+hankleSize+k-1,:); 
        end
        Yf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Yf;
        Uf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Uf;
        Up_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Up;
    end
    data_matrix = [Uf_tot Up_tot Yf_tot];
    [~ , R] = qr(data_matrix);
    L = R';
    L32 = L(2*(maxOrderNLE)*hankleSize+1:2*(maxOrderNLE)*hankleSize+hankleSize,(maxOrderNLE)*hankleSize+1:(maxOrderNLE)*hankleSize+maxOrderNLE*hankleSize);
    [Un,Sn,~] = svd(L32); 
    Sn = diag(Sn); 
    Sn = Sn(1:hankleSize); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hankleSize);
    if strcmp('preset', orderSelectMethodLE),
        m = orderLE;
    else
        m = orderselect(Sn,orderSelectMethodLE);
    end
    if m==0 
        condition = 0;
        warning('Selected order is zero.')
        warning('Attempt to identify the system failed.')
    end
    [AT , CT] = destac(R,m);
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Identified system is unstable.')
        warning('Attempt to identify the system failed.')
        condition = 0;
    end
else
    warning('Not enough number of samples is available for identification')
    warning('Attempt to identify the system failed.')
    condition = 0;
end
if condition>0
%Defining regressor matrices
%Gamma is the regressor for the initial conditions
    Gamma_total = zeros(size(out,1),p*m);
    Phi_total = zeros(size(out,1),(m+1)*maxOrderNLE);
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
        Phi = BD_omega_regressor(u(switch_time(i):switch_time(i+1)-1,:),AT,CT);
        Phi_total(switch_time(i):switch_time(i+1)-1,:) = Phi;
    end
    Phi = [Gamma_total Phi_total];
    %Iterative routine for static-nl, B, D and initial conditions estimation
    it=20;
    %b_hat is a vector whose first pm values are initial conditions and the
    %last m+1 values are B and D elements
    bd_hat = zeros(p * m+ m + 1,it);
    omega_hat = zeros(maxOrderNLE+1,it);
    omega0 = rand(maxOrderNLE+1,1);
    omega0 = omega0 / norm(omega0);
	omega0 = [0.01;halfwave_rectifier_tchebychev(min(in),max(in),maxOrderNLE-1)];
    omega0 = omega0 / norm(omega0);
    s1 = 10^10;
    s2 = 10^10;
    for i=1:it
        if i==1
            Phi_omega=[];
            for j = 2 : maxOrderNLE + 1
                temp_matrix = eye(m+1)*omega0(j);
                Phi_omega = [Phi_omega;temp_matrix];
            end
            temp_matrix = [eye(p * m)*omega0(1) zeros(p * m , m + 1)];
            temp_matrix_2 = zeros(maxOrderNLE * (m + 1), p * m);
            Phi_omega_final = [temp_matrix; temp_matrix_2 Phi_omega];
            Phi_omega_it = Phi*Phi_omega_final;
            bd_hat(:,i) = lscov(Phi_omega_it,out);
            sse_c = out'*out-bd_hat(:,i)'*Phi_omega_it'*out;
        else
            Phi_omega = [];
            for j = 2 : maxOrderNLE + 1
                temp_matrix = eye(m+1)*omega_hat(j,i - 1);
                Phi_omega = [Phi_omega;temp_matrix];
            end
            temp_matrix = [eye(p * m)*omega_hat(maxOrderNLE + 1 , i - 1) zeros(p * m , m + 1)];
            temp_matrix_2 = zeros(maxOrderNLE * (m + 1), p * m);
            Phi_omega_final = [temp_matrix; temp_matrix_2 Phi_omega];
            Phi_omega_it = Phi*Phi_omega_final;
            bd_hat(:,i) = lscov(Phi_omega_it,out);
            sse_c = out'*out-bd_hat(:,i)'*Phi_omega_it'*out;
        end
        Phi_bd = [];
        for j = 1 : maxOrderNLE
            temp_matrix = zeros(m+1,maxOrderNLE);
            temp_matrix(:,j) = bd_hat(p * m + 1:end,i);
            Phi_bd = [Phi_bd;temp_matrix];
        end
        temp_matrix = zeros(p * m , maxOrderNLE + 1);
        temp_matrix(:,1) = bd_hat(1 : p * m,i);
        temp_matrix_2 = zeros(maxOrderNLE * (m+1),1);
        Phi_bd_final = [temp_matrix;temp_matrix_2 Phi_bd];
        Phi_bd_it = Phi*Phi_bd_final;
        omega_hat(:,i) = lscov(Phi_bd_it,out);
        sse_b = out'*out-omega_hat(:,i)'*Phi_bd_it'*out;
        h = sign(omega_hat(1,i));
        bd_hat(:,i) = h*bd_hat(:,i)*norm(omega_hat(1:maxOrderNLE,i));
        omega_hat(:,i) = omega_hat(:,i)/norm(omega_hat(1:maxOrderNLE,i))*h;
        if (s1-sse_c<threshNSE) && (s2-sse_b<threshNSE)
            break
        end
            s1 = sse_c;
            s2 = sse_b;
    end
    it = i;
    %disp(['Terminated at iteration ',num2str(it)]);
    BT = bd_hat(p * m + 1:end-1,it);
    DT = bd_hat(end,it);
    DT = DT';
    initial = reshape(bd_hat(1 : p * m,it),m,length(bd_hat(1 : p * m,i))/m);
    initial = initial * omega_hat(1,it);
    omega = omega_hat(2:end,it);
    system_ss = ssm;
    set(system_ss,'A',AT,'B',BT,'C',CT,'D',DT,'domainIncr',ts,'nDelayInput',nDelayInput);
    newMin = min(in);
    newMax = max(in);
    newMean = mean(in);
    newStd = std(in);
    omega_coef = omega(:);
    static_nl = polynom('polyCoef',omega_coef,'polyType','tcheb','comment','Static Nonlinearity','polyRange',[newMin;newMax],'polyMean',newMean,'polyStd',newStd);
    system_nlbl = nlbl;
    set(system_nlbl,'elements',{static_nl,system_ss},'idMethod','subspace');
    if displayFlag == 1
        BT_kron = kron(BT,omega_coef');
        DT_kron = kron(DT,omega_coef');
        for i = 1 : p
            outp(switch_time(i):switch_time(i+1)-1) = dlsim(AT,BT_kron,CT,DT_kron,u(switch_time(i):switch_time(i+1)-1,:),initial(:,i));
        end
        outp = outp - mean(outp);
        if (plotMode == 1)
            h = figure;
            for i =1 : p
                figure(floor((i-1)/4)+h)
                subplot(4,1,mod(i-1,4)+1)            
                measured_data = out(switch_time(i):switch_time(i+1)-1);
                predicted_data = outp(switch_time(i):switch_time(i+1)-1);
                predicted_data = nldat(predicted_data','domainIncr',ts);
                measured_data = nldat(measured_data,'domainIncr',ts);
                set(measured_data,'chanNames','Measured output');
                set(predicted_data,'chanNames','Predicted output');
                plot(cat(2,measured_data,predicted_data),'plotmode','super');
            end
        end
    end
else
    static_nl = polynom;
    system_nlbl = nlbl;
    system_ss = ssm;
    set(system_nlbl,'elements',{static_nl,system_ss},'idMethod','subspace', ...
        'ordeLE',n);
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

function alpha = halfwave_rectifier_tchebychev(in_min,in_max,order)
x = in_min:0.0001:in_max;
y = max(x,0);
x = nldat(x','domainIncr',0.001);
y = nldat(y','domainIncr',0.001);
z = cat(2,x,y);
p = polynom(z,'polyType','tcheb','polyOrderMax',order,'polyOrderSelectMode','full');
alpha = p.polyCoef;
end
function dhat = intrinsicEstimator (g,k,y)
    Hg=(eye(size(g,2))-pinv(g)*k*pinv(k)*g);
    Gg=pinv(g)-pinv(g)*k*pinv(k);
    dhat = pinv(Hg)*Gg*y;
end
