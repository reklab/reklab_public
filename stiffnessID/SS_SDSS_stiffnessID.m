function system = SS_SDSS_stiffnessID (z,varargin)
% system= pcas_short_segment (z,switch_time,varargin)
% This function requires NLID toolbox in the MATLAB path
% This function estimates parallel-cascade joint stiffness between input and output
%
%
%
options={{'decimation_ratio' 10 'decimation ratio'} ...
         {'maxordernle' 8 'maximum order for nonlinearity'} ...
         {'hanklesize' 15 'Size of hankle matrix'} ...
         {'delayinput' 0.05 'Delay added to the input'} ...
         {'orderselectmethod','manual'}...
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
% global tqI_all_segments tqR_all_segments
% tqi = tqI_all_segments.dataSet;
% tqi = tqi - mean(tqi);
% tqr = tqR_all_segments.dataSet;
% tqr = tqr - mean(tqr);
% tqI_noisy = zeros(sum(segLength),1);
% tqR_noisy = zeros(sum(segLength),1);
pointer = 1;
switch_time = zeros(length(endpointer)-1,1);
for i = 1 : length(endpointer)
    vel_seg = velocity(onsetPointer(i):endpointer(i));
    dvel_seg = del(vel_seg,delayinput/ts);
    dvelocity(pointer:pointer+segLength(i)-1) = dvel_seg;
    tqT_noisy(pointer:pointer+segLength(i)-1) = output(onsetPointer(i):endpointer(i));
%     tqI_noisy(pointer:pointer+segLength(i)-1) = tqi(onsetPointer(i):endpointer(i));
%     tqR_noisy(pointer:pointer+segLength(i)-1) = tqr(onsetPointer(i):endpointer(i));
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
% [tqI_noisy,~,~,~] = decimate_segment(tqI_noisy,switch_time(1:end-1),decimation_ratio);
% [tqR_noisy,~,~,~] = decimate_segment(tqR_noisy,switch_time(1:end-1),decimation_ratio);
[tqT_noisy,switch_time,segLength,~] = decimate_segment(tqT_noisy,switch_time(1:end-1),decimation_ratio);
% tqI_noisy = tqI_noisy - mean(tqI_noisy);
% tqR_noisy = tqR_noisy - mean(tqR_noisy);
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
    [Un,Sn,~] = svd(L32); 
    Sn = diag(Sn); 
    Sn = Sn(1:hanklesize); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hanklesize);
%     m = orderselect(Sn,orderselectmethod);
    m = 2;
    [AT , CT] = destac(R,m);
    if m==0 
        condition = 0;
        warning('Selected reflex system order is set to zero.')
        warning('Only the intrinsic pathway will be identified.')
    end
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Reflex linear system is unstable.')
        warning('Attempt to identify a reflex pathway failed.')
        warning('Only the intrinsic pathway will be identified.')
        condition = 0;
    end
else
    warning('Not enough number of samples is available')
    warning('Attempt to identify a reflex pathway failed.')
    warning('Only the intrinsic pathway will be identified.')
    condition = 0;
end
%Identify the intrinsic pathway independently than estimate of reflexes
if condition>0
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
    intrinsic = opt_li_wen_2011 (u_i,Phi,tqT_noisy);
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
    [Un,Sn,~] = svd(L32); 
    Sn = diag(Sn); 
    Sn = Sn(1:hanklesize); 
    R = struct('L',L,'Un',Un,'m',1,'l',1,'i',hanklesize);
    [AT , CT] = destac(R,m);
    if ~isempty(find(abs(eig(AT))>1, 1))
        warning('Reflex linear system is unstable.')
        warning('Attempt to identify a reflex pathway failed.')
        warning('Only the intrinsic pathway will be identified.')
        condition = 0;
    end
    if condition > 0
%Iterative routine for static-nl, B, D and initial conditions estimation
        it=20;
%b_hat is a vector whose first pm values are initial conditions and the
%last m+1 values are B and D elements
        bd_hat = zeros(p * m+ m + 1,it);
        omega_hat = zeros(maxordernle+1,it);
        %omega0 = ones(order+1,1);
        omega0 = [0.01;halfwave_rectifier_tchebychev(min(dvelocity),max(dvelocity),maxordernle-1)];
        %omega0 = halfwave_rectifier_tchebychev(-1,+1,order);
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
        v = vaf(tqT_noisy-tqI,tqR);
        v = v.dataSet;
        %disp(['VAF between intrinsic residual and reflex torque: ',num2str(v)])
%        tqI_noisy = nldat(tqI_noisy,'domainIncr',ts);
%        tqR_noisy = nldat(tqR_noisy,'domainIncr',ts);
%         figure
%         subplot(4,2,1)
%         i=2;
%         predicted_data = u_i(switch_time(i):switch_time(i+1)-1,6);
%         predicted_data = predicted_data;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         set(predicted_data,'chanNames','Position (rad)');
%         plot(predicted_data);
%         title('(A)')
%         subplot(4,2,2)
%         i=4;
%         predicted_data = u_i(switch_time(i):switch_time(i+1)-1,6);
%         predicted_data = predicted_data;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         set(predicted_data,'chanNames','Position (rad)');
%         plot(predicted_data);
%         title('(B)')
%         subplot(4,2,3)
%         i=2;
%         measured_data = tqT_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqT(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         measured_data = nldat(measured_data,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(C)')
%         subplot(4,2,4)
%         i=4;
%         measured_data = tqT_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqT(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         measured_data = nldat(measured_data,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(D)')
%         subplot(4,2,5)
%         i=2;
%         measured_data = tqI_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqI(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         measured_data = nldat(measured_data,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(E)')
%         subplot(4,2,6)
%         i=4;
%         measured_data = tqI_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqI(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data,'domainIncr',ts);
%         measured_data = nldat(measured_data,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(F)')
%         subplot(4,2,7)
%         i=2;
%         measured_data = tqR_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqR(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data+4.3,'domainIncr',ts);
%         measured_data = nldat(measured_data+4.3,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(G)')
%         subplot(4,2,8)
%         i=6;
%         measured_data = tqR_noisy(switch_time(i):switch_time(i+1)-1);
%         measured_data = measured_data.dataSet;
%         measured_data = measured_data - mean(measured_data);
%         predicted_data = tqR(switch_time(i):switch_time(i+1)-1);
%         predicted_data = predicted_data.dataSet;
%         predicted_data = predicted_data - mean(predicted_data);
%         predicted_data = nldat(predicted_data+4.3,'domainIncr',ts);
%         measured_data = nldat(measured_data+4.3,'domainIncr',ts);
%         set(measured_data,'chanNames','Measured torque');
%         set(predicted_data,'chanNames','Predicted torque');
%         plot(cat(2,measured_data,predicted_data),'plotmode','super');
%         hold on
%         plot(measured_data-predicted_data,'line_color','r')
%         title('(H)')
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
    condition = 0;
end
if (condition==0)
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
        pos_seg = input(onsetPointer(i):endpointer(i));
        pos_seg = nldat(pos_seg,'domainIncr',ts);
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