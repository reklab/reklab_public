function system = SS_SDSS_stiffnessID (z,varargin)
% system = SS_SDSS_stiffnessID (z)
% This function requires NLID toolbox in MATLAB path
% This function estimates parallel-cascade joint stiffness between input and output
% from short segments of data
%
%
options={{'decimation_ratio' 10 'decimation ratio'} ...
         {'maxordernle' 8 'maximum order for nonlinearity'} ...
         {'hanklesize' 20 'Size of hankle matrix'} ...
         {'delayinput' 0.05 'Delay added to the input'} ...
         {'orderselectmethod','auto'}...
         {'threshold' 10^(-5)}...
         {'plot_mode' 0 '1 to plot and 0 to not plot segments'}...
     };
if arg_parse(options,varargin);
     return
 end
% Author: Kian Jalaleddini
% Date: February 11, 2013 Ver 0.1
% Date: September 17, 2013 Ver 0.2
% Date: October 9, 2013 Ver 0.3 % adding refinement steps of the A & C matrices
% Date: March 1, 2014 Ver 0.4 adding compatibility with segdat objects
% Date: April 29, 2014 Ver 0.5 adding irf for the intrinsic pathway
% Date: May 12, 2014 Ver. 0.6 correcting issues with input initial conditions
%%
reflexPathID = 1;
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
input = input - mean(input);
output = data(:,2);
output = output - mean(output);
N = floor(segLength/decimation_ratio) - 2 * hanklesize + 1;
if length(find(N<1))>1
    warning(['Removing ',num2str(length(find(N<1))),' very short segments'])
end
onsetPointer(N<1) = [];
segLength(N<1) = [];
N(N<1) = [];
if ~isempty(N)
endpointer = onsetPointer + segLength - 1;
%extracting input-output data from segdat
irf_len_i = delayinput/ts/decimation_ratio;
lags_i = (-irf_len_i:1:irf_len_i);
nLags_i = length(lags_i);
%inputDec = decimate(input,decimation_ratio);
positionDelay = zeros(size(input,1),nLags_i);
for j = 1:nLags_i
    positionDelay(:,j) = del(input,lags_i(j)*decimation_ratio);
end
velocity = ddt(nldat(input,'domainIncr',ts));
velocity = get(velocity,'dataSet');
positionDelaysegments = zeros(sum(segLength),nLags_i);
dvelocity = zeros(sum(segLength),1);
tqT_noisy = zeros(sum(segLength),1);
pointer = 1;
switch_time = zeros(length(endpointer)-1,1);
for i = 1 : length(endpointer)
    vel_seg = velocity(onsetPointer(i):endpointer(i));
    dvel_seg = del(vel_seg,delayinput/ts);
    dvelocity(pointer:pointer+segLength(i)-1) = dvel_seg;
    tqT_noisy(pointer:pointer+segLength(i)-1) = output(onsetPointer(i):endpointer(i));
    positionDelaysegments(pointer:pointer+segLength(i)-1,:) = positionDelay(onsetPointer(i):endpointer(i),:);
    pointer = pointer + segLength(i);
    switch_time(i) = pointer;
end
[dvelocity,~,~,~] = decimate_segment(dvelocity,switch_time(1:end-1),decimation_ratio);
positionDelaysegments = bsxfun(@minus,positionDelaysegments,mean(positionDelaysegments));
u_i = zeros(size(dvelocity,1),nLags_i);
for i = 1:nLags_i
    [u_i(:,i),~,~,~] = decimate_segment(positionDelaysegments(:,i),switch_time(1:end-1),decimation_ratio);
end
[tqT_noisy,switch_time,segLength,~] = decimate_segment(tqT_noisy,switch_time(1:end-1),decimation_ratio);
tqT_noisy = tqT_noisy - mean(tqT_noisy);

N = segLength - 2 * hanklesize + 1;
ts = ts * decimation_ratio;
nsamp = length(dvelocity);
p = length(segLength);
%Ensure enough number of samples is available O.W. identify intrinsic path only
if nsamp>2*hanklesize*p-p+2*maxordernle*hanklesize+nLags_i*hanklesize+1
%Construct the input signal
    avg = (max(dvelocity) + min(dvelocity)) / 2;
    rng = max(dvelocity) - min(dvelocity);
    un = (dvelocity - avg) * 2 / rng;
    u_r = multi_tcheb(un,maxordernle - 1);
    u = [u_i,u_r];
%First attempt to identify AT and CT
    Yf_tot = zeros(sum(N),hanklesize);
    Uf_tot = zeros(sum(N),(maxordernle+nLags_i)*hanklesize);
    Up_tot = zeros(sum(N),(maxordernle+nLags_i)*hanklesize);
    for i = 1 : p
        Uf = zeros(N(i),(maxordernle+nLags_i)*hanklesize); 
        Up = zeros(N(i),(maxordernle+nLags_i)*hanklesize); 
        Yf = zeros(N(i), hanklesize);
        u_r_segment = u(switch_time(i):switch_time(i+1)-1,:);
        output_segment = tqT_noisy(switch_time(i):switch_time(i+1)-1,:);
        for k = (1:hanklesize)
            Up(:,(k-1) * (maxordernle+nLags_i)+1:k * (maxordernle+nLags_i)) = u_r_segment(k:N(i)+k-1,:); 
            Uf(:,(k-1) * (maxordernle+nLags_i)+1:k * (maxordernle+nLags_i)) = u_r_segment(hanklesize+k:N(i)+hanklesize+k-1,:); 
            Yf(:,(k-1) * 1+1:k * 1) = output_segment(hanklesize+k:N(i)+hanklesize+k-1,:); 
        end
        Yf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Yf;
        Uf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Uf;
        Up_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Up;
    end
    data_matrix = [Uf_tot Up_tot Yf_tot];
    [~ , R] = qr(data_matrix);
    L = R';
    L32 = L(2*(maxordernle+nLags_i)*hanklesize+1:2*(maxordernle+nLags_i)*hanklesize+hanklesize,(maxordernle+nLags_i)*hanklesize+1:(maxordernle+nLags_i)*hanklesize+maxordernle*hanklesize);
    [Un,~,~] = svd(L32); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hanklesize);
    %m = orderselect(Sn,orderselectmethod);
    m = 2;
    [AT , CT] = destac(R,m);
    if m==0 
        reflexPathID = 0;
        warning('Selected reflex system order is set to zero.')
        warning('Only the intrinsic pathway will be identified.')
    end
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Reflex linear system is unstable.')
        warning('Attempt to identify a reflex pathway failed.')
        warning('Only the intrinsic pathway will be identified.')
        reflexPathID = 0;
    end
else
    warning('Not enough number of samples is available')
    warning('Attempt to identify a reflex pathway failed.')
    warning('Only the intrinsic pathway will be identified.')
    reflexPathID = 0;
end
%Identify the intrinsic pathway independently than estimate of reflexes
if reflexPathID>0
%Defining regressor matrices
%Gamma is the regressor for the initial conditions
    Gamma_total = zeros(size(tqT_noisy,1),p*m);
    Phi_total = zeros(size(tqT_noisy,1),(m+1)*maxordernle);
    max_interval = max(segLength);
    Gamma_nominal = zeros(max_interval,m);
    Gamma_nominal(1,:) = CT;
    An = AT;
    for i = 1:floor(log(nsamp)/log(2))
        Gamma_nominal(2^(i-1)+1:2^i,:) = Gamma_nominal(1:2^(i-1),:)*An;
        An = An * An;
    end
    Gamma_nominal(2^i+1:nsamp,:) = Gamma_nominal(1:nsamp-2^i,:) * An;
    for i = 1 : p
        Gamma_total(switch_time(i):switch_time(i+1)-1,(i-1)*m+1:i*m) = Gamma_nominal(1:segLength(i),:);
        Phi = BD_omega_regressor(u_r(switch_time(i):switch_time(i+1)-1,:),AT,CT);
        Phi_total(switch_time(i):switch_time(i+1)-1,:) = Phi;
    end
    Phi = [Gamma_total Phi_total];
    intrinsic = intrinsicEstimator(u_i,Phi,tqT_noisy);
    tqI = u_i * intrinsic;
    tqI_res = tqT_noisy - tqI;
    tqI_res = tqI_res - mean(tqI_res);

%Second attempt to refine the estimates of A and C
    Yf_tot = zeros(sum(N),hanklesize);
    Uf_tot = zeros(sum(N),(maxordernle)*hanklesize);
    Up_tot = zeros(sum(N),(maxordernle)*hanklesize);
    for i=1:p
        Uf = zeros(N(i),(maxordernle) * hanklesize); 
        Up = zeros(N(i),(maxordernle) * hanklesize); 
        Yf = zeros(N(i), hanklesize);
        u_r_segment = u_r(switch_time(i):switch_time(i+1)-1,:);
        output_segment = tqI_res(switch_time(i):switch_time(i+1)-1,:);
        for k = (1:hanklesize)
          Up(:,(k-1) * maxordernle+1:k * maxordernle) = u_r_segment(k:N(i)+k-1,:); 
          Uf(:,(k-1) * maxordernle+1:k * maxordernle) = u_r_segment(hanklesize+k:N(i)+hanklesize+k-1,:); 
          Yf(:,(k-1) * 1+1:k * 1) = output_segment(hanklesize+k:N(i)+hanklesize+k-1,:); 
        end
        Yf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Yf;
        Uf_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Uf;
        Up_tot(sum(N(1:i))-N(i)+1:sum(N(1:i)),:) = Up;
    end
    data_matrix = [Uf_tot Up_tot Yf_tot];
    [~ , R] = qr(data_matrix);
    L = R';
    L32 = L(2*maxordernle*hanklesize+1:2*maxordernle*hanklesize+hanklesize,maxordernle*hanklesize+1:2*maxordernle*hanklesize);
    [Un,~,~] = svd(L32); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hanklesize);
    [AT , CT] = destac(R,m);
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Reflex linear system is unstable.')
        warning('Attempt to identify a reflex pathway failed.')
        warning('Only the intrinsic pathway will be identified.')
        reflexPathID = 0;
    end
    if reflexPathID > 0
%Iterative routine for static-nl, B, D and initial conditions estimation
        it=20;
%b_hat is a vector whose first pm values are initial conditions and the
%last m+1 values are B and D elements
        bd_hat = zeros(p * m+ m + 1,it);
        omega_hat = zeros(maxordernle+1,it);
        omega0 = [0.01;halfwave_rectifier_tchebychev(min(dvelocity),max(dvelocity),maxordernle-1)];
        omega0 = omega0 / norm(omega0);
        s1 = 10^10;
        s2 = 10^10;
        for i=1:it
            if i==1
                Phi_omega=[];
                for j = 2 : maxordernle + 1
                    temp_matrix = eye(m+1)*omega0(j);
                    Phi_omega = [Phi_omega;temp_matrix];
                end
                temp_matrix = [eye(p * m)*omega0(1) zeros(p * m , m + 1)];
                temp_matrix_2 = zeros(maxordernle * (m + 1), p * m);
                Phi_omega_final = [temp_matrix; temp_matrix_2 Phi_omega];
                Phi_omega_it = Phi*Phi_omega_final;
                bd_hat(:,i) = lscov(Phi_omega_it,tqI_res);
                sse_c = tqI_res'*tqI_res-bd_hat(:,i)'*Phi_omega_it'*tqI_res;
            else
                Phi_omega = [];
                for j = 2 : maxordernle + 1
                    temp_matrix = eye(m+1)*omega_hat(j,i - 1);
                    Phi_omega = [Phi_omega;temp_matrix];
                end
                temp_matrix = [eye(p * m)*omega_hat(maxordernle + 1 , i - 1) zeros(p * m , m + 1)];
                temp_matrix_2 = zeros(maxordernle * (m + 1), p * m);
                Phi_omega_final = [temp_matrix; temp_matrix_2 Phi_omega];
                Phi_omega_it = Phi*Phi_omega_final;
                bd_hat(:,i) = lscov(Phi_omega_it,tqI_res);
                sse_c = tqI_res'*tqI_res-bd_hat(:,i)'*Phi_omega_it'*tqI_res;
            end
            Phi_bd = [];
            for j = 1 : maxordernle
                temp_matrix = zeros(m+1,maxordernle);
                temp_matrix(:,j) = bd_hat(p * m + 1:end,i);
                Phi_bd = [Phi_bd;temp_matrix];
            end
            temp_matrix = zeros(p * m , maxordernle + 1);
            temp_matrix(:,1) = bd_hat(1 : p * m,i);
            temp_matrix_2 = zeros(maxordernle * (m+1),1);
            Phi_bd_final = [temp_matrix;temp_matrix_2 Phi_bd];
            Phi_bd_it = Phi*Phi_bd_final;
            omega_hat(:,i) = lscov(Phi_bd_it,tqI_res);
            sse_b = tqI_res'*tqI_res-omega_hat(:,i)'*Phi_bd_it'*tqI_res;
            h = sign(omega_hat(1,i));
            bd_hat(:,i) = h*bd_hat(:,i)*norm(omega_hat(1:maxordernle,i));
            omega_hat(:,i) = omega_hat(:,i)/norm(omega_hat(1:maxordernle,i))*h;
            if (s1-sse_c<threshold) && (s2-sse_b<threshold)
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
        system_ss = ss(AT,BT,CT,DT,ts);
        tf_l = tf(system_ss);
        num = get(tf_l,'num');
        num = num{1};
        den = get(tf_l,'den');
        den = den{1};
        gain = sum(num)/sum(den);
        system_ss = ssm;
        set(system_ss,'A',AT,'B',BT/gain,'C',CT,'D',DT/gain,'domainIncr',ts,'nDelayInput',delayinput/ts);
        newMin = min(dvelocity);
        newMax = max(dvelocity);
        newMean = mean(dvelocity);
        newStd = std(dvelocity);
        omega_coef = omega(:);
        static_nl = polynom('polyCoef',omega*gain,'polyType','Tcheb','comment','Static Nonlinearity','polyRange',[newMin;newMax],'polyMean',newMean,'polyStd',newStd);
        tqR = zeros(size(tqT_noisy));
        BT_kron = kron(BT,omega_coef');
        DT_kron = kron(DT,omega_coef');
        for i = 1 : p
            tqR(switch_time(i):switch_time(i+1)-1) = dlsim(AT,BT_kron,CT,DT_kron,u_r(switch_time(i):switch_time(i+1)-1,:),initial(:,i));
        end
        tqR = nldat(tqR,'domainIncr',ts);
        tqI = nldat(tqI,'domainIncr',ts);
        tqT = tqI + tqR;
        tqT_noisy = nldat(tqT_noisy,'domainIncr',ts);
        vaf_tot = vaf(tqT_noisy,tqT);
        vaf_I = vaf(tqT_noisy,tqI);
        vaf_R = vaf(tqT_noisy,tqR);
        if plot_mode == 1
            for i =1 : p
                figure(floor((i-1)/4)+10)
                subplot(4,1,mod(i-1,4)+1)
                measured_data = tqT_noisy(switch_time(i):switch_time(i+1)-1);
                measured_data = measured_data.dataSet;
                measured_data = measured_data - mean(measured_data);
                predicted_data = tqT(switch_time(i):switch_time(i+1)-1);
                predicted_data = predicted_data.dataSet;
                predicted_data = predicted_data - mean(predicted_data);
                predicted_data = nldat(predicted_data,'domainIncr',ts);
                measured_data = nldat(measured_data,'domainIncr',ts);
                set(measured_data,'chanNames','Measured torque');
                set(predicted_data,'chanNames','Predicted torque');
                plot(cat(2,measured_data,predicted_data),'plotmode','super');
                hold on
                plot(measured_data-predicted_data,'line_color','r')
            end
        end
    end
end
else
    reflexPathID = 0;
end
if (reflexPathID==0)
%Attempt to estimate the reflex path failed, only estimate the intrinsic path
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
    %input = input - mean(input);
    output = data(:,2);
    output = output - mean(output);
    endpointer = onsetPointer + segLength - 1;
    %extracting input-output data from segdat
    positionDelay = zeros(size(input),nLags_i);
    for j = 1:nLags_i
        positionDelay(:,j) = del(input,lags_i(j)*decimation_ratio);
    end
    positionDelaysegments = zeros(sum(segLength),nLags_i);
    tqT_noisy = zeros(sum(segLength),1);
    pointer = 1;
    switch_time = zeros(length(endpointer)-1,1);
    irf_len_i = delayinput/ts/decimation_ratio;
    lags_i = (-irf_len_i:1:irf_len_i);
    nLags_i = length(lags_i);
    for i = 1 : length(endpointer)
        positionDelaysegments(pointer:pointer+segLength(i)-1,:) = positionDelay(onsetPointer(i):endpointer(i),:);
        tqT_noisy(pointer:pointer+segLength(i)-1) =output(onsetPointer(i):endpointer(i));
        pointer = pointer + segLength(i);
        switch_time(i) = pointer;
    end
    positionDelaysegments = bsxfun(@minus,positionDelaysegments,mean(positionDelaysegments));
    [tqT_noisy1,~,~,~] = decimate_segment(tqT_noisy,switch_time(1:end-1),decimation_ratio);
    u_i = zeros(size(tqT_noisy1,1),nLags_i);
    for i = 1:nLags_i
        [u_i(:,i),~,~,~] = decimate_segment(positionDelaysegments(:,i),switch_time(1:end-1),decimation_ratio);
    end
    [tqT_noisy,~,~,~] = decimate_segment(tqT_noisy,switch_time(1:end-1),decimation_ratio);
    tqT_noisy = tqT_noisy - mean(tqT_noisy);
    ts = ts * decimation_ratio;
    intrinsic=u_i\tqT_noisy;
    tqI = nldat(u_i*intrinsic,'domainIncr',ts);
    tqR = tqI*0;
    tqT = tqI;
    tqT_noisy = nldat(tqT_noisy,'domainIncr',ts);
    tqI = nldat(tqI,'domainIncr',ts);
    %tqR = nldat(tqR,'domainIncr',ts);
    vaf_tot = vaf(tqT_noisy,tqT);
    vaf_I = vaf(tqT_noisy,tqI);
    vaf_R = vaf(tqT_noisy,tqR);
    system_ss = ssm;
    set(system_ss,'domainIncr',ts,'nDelayInput',delayinput/ts);
    static_nl = polynom;
end
%Assigning Function's output
reflex = cell(2,1);
reflex{1} = static_nl;
reflex{2} = system_ss;
vafs = [vaf_tot.dataSet;vaf_I.dataSet;vaf_R.dataSet];
vafs((vafs>100)) = 0;
vafs((vafs<0)) = 0;
system = cell(3,1);
intrinsic = irf('nSides',2,'dataSet',intrinsic/ts,'domainIncr',ts,'domainStart',-irf_len_i*ts,'comment','Intrinsic IRF','chanNames','IRF');
system{1} = intrinsic;
system{2} = reflex;
system{3} = vafs;
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
function [output,switch_time_new,interval,p] = decimate_segment(input,switch_time,decimation_ratio)
sw = [1;switch_time];
p = length(sw);
sw = [sw; length(input)+1];
output = [];
switch_time_new = 1;
interval=zeros(length(switch_time),1);
segment_onset=1;
for i = 1 : p
    output_temp = decimate(input(sw(i):sw(i+1)-1),decimation_ratio);
    output = [output;output_temp];
    segment_onset = segment_onset+length(output_temp);
    switch_time_new = [switch_time_new;segment_onset];
    interval(i) = length(output_temp);
end
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
function [d_x] = del(x,nDelay)
% this function adds a delay to the input signal 
% x is the input signal
% nDelay - dealy in samples
% d_x is the delayed singal
% Now supports negative delays
d_x = zeros(size(x));
if nDelay>=0
    d_x(nDelay+1:end) = x(1:end-nDelay);
elseif nDelay<0
    nDelay = abs(nDelay);
    d_x(1:end-nDelay) = x(nDelay+1:end);
end
end
function [W, flags] = multi_tcheb(V, max_order);
%
%  usage W = multi_herm2(V, max_order);
%
%  given a collection of column vectors, V, this function returns
%  a matrix of all of the Tcheb functions up to max_order applied
%  to all of the vectors in V.

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

if max_order==0,
   W=V*0 +1;
   return;
end




[nr,nc] = size(V);


%  create a matrix of flags, such that flags (i,j) points to the 
%  column of the i'th basis function raised to the j'th power.

flags = zeros(nc,max_order);
flags(:,1) = [2:nc+1]';
for order = 2:max_order
   % first column is offsetr by 1 from last value of previous order
   flags(1,order) = flags(nc,order-1)+1;
    for i = 2 : nc
        num_terms = flags(nc,order-1)-flags(i-1,order-1)+1;
        flags(i,order) = flags(i-1,order)+ num_terms;
      end
   end

   

% pre-allocate the W matrix
W = zeros(nr,flags(nc,max_order));
%  generate the single basis function Hermite polynomials and place 
%  them  in the correct columns of W

Te = zeros(nr,max_order+1);
Te(:,1) = ones(nr,1);

%  generate all of the functions involving a Hermite polynomial
%  applied to a single basis vector

for i = 1:nc
    Te(:,2) = V(:,i);
    for j = 2:max_order
        Te(:,j+1) = 2*V(:,i).*Te(:,j) - Te(:,j-1);
      end
    for j = 1:max_order
        W(:,flags(i,j)) = Te(:,j+1);
      end
  end
W(:,1) = ones(nr,1);

clear Te V


%  generate all of the functions involving a powers of a single input vector


  
%  Now, using the functions that we just created, fill in the 
%  rest of the matrix



for order = 2:max_order
    for v1 = 1 : nc-1
        index = flags(v1,order);
        for v1_order = order-1:-1:1
            term1 = W(:,flags(v1,v1_order));
            rem_order = order - v1_order;
            
%           find the terms of order rem_order, whose leading variable
%           is 'greater' than v1

            first_term = flags(v1+1,rem_order);
            last_term = flags(nc,rem_order);
            for j = first_term:last_term
                index = index+1;
                W(:,index)=term1.*W(:,j);
                disp([index flags(v1,v1_order) j])
              end
          end
      end
  end
end



