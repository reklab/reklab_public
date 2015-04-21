function [intrinsic, reflex, tqI, tqR, tqT] = SDSS_stiffnessID (z,varargin)
% [system_ss, system, K, B, I, tqR, tqI, tqT, error]= pcascsub (z)
% This function estimates a parallel-cascade model between input and output stored in columns of z.
% z = cat(2 , position, torque)
% The output, pc_nss is a 3x1 cell
% pc{1} = [K;B;I];
% pc{2} = reflex; where reflex{1}=static nonlinearity, reflex{2}=state-space model
% pc{3} = vafs; where vafs=[v]
% This function is based on the following paper:
%[*] K. Jalaleddini, Ehsan Sobhani Tehrani and R. E. Kearney, "A Subspace Approach to the Structural Decomposition and Identification of Ankle Joint Dynamic Stiffness", IEEE TBME.
%[1] K. Jalaleddini and R. E. Kearney, " Subspace Identification of SISO Hammerstein Systems: Application to Stretch Reflex Identification", IEEE TBME 2013.
options={{'decimation_ratio' 10 'ratio to decimate data'} ...
         {'order' 12 'maximum order for nonlinearity'} ...
         {'hankle_size' 20 'Size of hankle matrix'} ...
         {'delay' 0.04 'Reflex delay in s'} ...
         {'orderdetection','manual'}...
     };
 if arg_parse(options,varargin);
     return
 end
% Author: Kian Jalaleddini 
%% Initialization
pos = z(:,1);
trq = z(:,2);
z = cat(2,pos,trq);
ts = get(z,'domainIncr');
pos = pos - mean(pos);
trq = trq - mean(trq);
%Define position, velocity and acceleration signals.
vel = ddt(pos);
pos = get(pos,'dataSet');
velData = get(vel,'dataSet');
vel = decimate(vel,decimation_ratio);
trq = get(trq,'dataSet');
dvel = del(velData,ceil(delay/(ts)));
dvel = decimate(dvel,decimation_ratio);
pos = decimate(pos,decimation_ratio);
trq = decimate(trq,decimation_ratio);
nSamp = size(pos,1);
%Constructed input signal based on (5) in [*] using Tchebycehv polynomial
avg = (max(dvel) + min(dvel))/2;
rng = max(dvel) - min(dvel);
un = (dvel - avg)*2/rng;
u_r = multi_tcheb(un,order-1);
%Constructed intrinsic regressor
irf_len_i = delay/ts/decimation_ratio;
lags_i = (-irf_len_i:1:irf_len_i);
nLags_i = length(lags_i);
u_i = zeros(length(pos),nLags_i);
for i = 1:nLags_i
    u_i(:,i) = del(pos,lags_i(i));
end
u = [u_i,u_r];
%% Identification
%Estimating AT and CT from SMI toolbox
[S, R] = dordpi(u,trq,hankle_size);
n = orderselect(S,orderdetection);
[AT, CT] = destac(R,n);
f=find(abs(eig(AT))>1, 1);
if ~isempty(f)
    warning('Attempt to identify a stable reflex pathway failed.')
    warning('Identified reflex linear system is unstable..')
    warning('Only intrinsic pathway will be identified.')
    n=0;
end
%If the order of the linear system is selected to be 0, then it only
%estimates the intrinsic pathway using IRF estimation.
if n>0
    %If the order of the reflex system is larger than 0
    %
    %
    %Estimating the regressor Phi based on (9) in [*]
    temp_matrix = zeros(nSamp,n*order);
    e=eye(n);
    for j = 1 : order
        for i=1:n
            x=ltitr(AT,e(:,i),u_r(:,j));
            yij=CT*x';
            temp_matrix(:,(j-1)*n+i)=yij(:);
        end
    end
    Phi = zeros(nSamp,(n+1)*order);
    k1 = 1;
    k2 = 1;
    for i = 1 : (n+1)*order
        if mod(i,n+1)==0
            Phi(:,i) = u_r(:,k1);
            k1 = k1+1;
        else
            Phi(:,i) = temp_matrix(:,k2);
            k2 = k2+1;
        end
    end
    Phi(1,:)=1;
    %Orthogonal projection to estimate intrinsic path (13) in [*]
    Hg = eye(length(lags_i)) - pinv(u_i) * Phi * pinv(Phi)* u_i;
    Gg = pinv(u_i) - pinv(u_i) * Phi * pinv(Phi);
    intrinsic = pinv(Hg) * Gg * trq;
    tqI = u_i * intrinsic;
    intrinsic = irf('nSides',2,'dataSet',intrinsic/(ts*decimation_ratio),'domainIncr',ts*decimation_ratio,'domainStart',-irf_len_i*ts*decimation_ratio,'comment','Identified intrinsic IRF','chanNames','IRF');
    tqI_r = trq - tqI;
    tqI = nldat(tqI,'domainIncr',ts*decimation_ratio);
    zReflex = cat(2,vel,nldat(tqI_r,'domainIncr',ts*decimation_ratio));
    reflex = nlbl(zReflex,'idMethod','subspace','nDelayInput',delay/ts/decimation_ratio,'maxOrderNLE',order,'threshNSE',10^-5,'displayFlag',0,'hankleSize',hankle_size,'orderSelect',orderdetection);
    set(reflex,'comment','Identified reflex Hammerstein');
    if isempty(reflex{2}.A)
        n = 0;
        warning('Attempt to identify a stable reflex pathway failed.')
        warning('Identified reflex linear system is unstable..')
        warning('Only intrinsic pathway will be identified.')
    else
        tqR = nlsim(reflex,vel);
        tqT = tqI + tqR;
    end
if n==0
    % If the order of the reflex system is selected as 0 only identify the
    % intrinsic pathway
    lags_i = (-irf_len_i:1:irf_len_i);
    nLags_i = length(lags_i);
    u_i = zeros(length(pos),nLags_i);
    for i = 1:nLags_i
        u_i(:,i) = del(pos,lags_i(i));
    end
    intrinsic=u_i\trq;
    tqI=nldat(u_i*intrinsic,'domainIncr',ts/decimation_ratio);
    intrinsic = irf('nSides',2,'dataSet',intrinsic/(ts*decimation_ratio),'domainIncr',ts*decimation_ratio,'domainStart',-irf_len_i*ts*decimation_ratio,'comment','Intrinsic IRF','chanNames','IRF');
    tqR=tqI*0;
    tqT=tqI;
    nl_reflex = polynom;
    l_reflex_ss =ssm;
    reflex = nlbl('elements',{nl_reflex,l_reflex_ss},'comment','Relex stiffness model');
end
end