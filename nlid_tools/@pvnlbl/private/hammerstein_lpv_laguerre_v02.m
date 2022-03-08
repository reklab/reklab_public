function [static_nl,h_r,yp,fIter]= hammerstein_lpv_laguerre_v02(z,rho,wNL,wIRF,varargin)
% This function estimates an LPV Hammerstein Model from Input/Output data.
% The model structure identified is an LPV Hammerstein (LPV Static NL + LPV Laguerre Expansion of IRF).
% The inputs (in the NLDAT object z) are:
% -- 1) z: Input/Output data
% -- 2) rho: The scheduling variable
% -- 3) wNL: Weights related to NL parameters: 1 if a parameter is required and 0 if a parameter must be excluded (added on Aug. 28, 2015 to accomodate NL parameter selection)
% -- 4) wIRF: Weights related to NL parameters: 1 if a parameter is required and 0 if a parameter must be excluded (added on Aug. 28, 2015 to accomodate Laguerre IRF parameter selection)

% The outputs are:
%--- 1) static_nl is the 2D static nonlinearity with dimension [n*p,1].
%       n is the order of expansion of delayed velocity and p is the order of expansion of the scheduling variable.

%--- 2) h_r: The structure for LPV Laguerre exapnsion of the reflex pathway IRF. The structure has the following parameters: 1)h_r.lpv_laguerre; 2)h_r.domainIncr
%       The dimension of h_r.data is [m*q,1], where q is the order of the Laguerre expansion of the reflex IRF and m is the expansion order of the reflex IRF Laguerre w.r.t. scheduling variable.

%--- 3) yp is the output predicted by the model
%--- 4) Added final iteration, fIter, at which the algorithm stopped/converged

options={{'decimation_ratio' 10 'ratio to decimate data bt'} ...
         {'irf_len_r' 0.6 'length of reflex IRF in s'} ...
         {'nside_r' 1 'number of sides of intrinsic IRF'} ...
         {'q' 14 'order of the Laguerre expansion of the reflex IRF'} ...
         {'alfa' 0.8 'constant parameter of the Laguerre expansion'} ...
         {'m' 2 'order of SV expansion for reflex IRF Laguerre exapnsion'} ...
         {'n' 8 'order of velocity expansion for static nonlinearity'} ...
         {'p' 3 'order of SV expansion for static nonlinearity'} ...
         {'max_iter' 500 'maximum number of iterations for reflex pathway identification'} ...
         {'nl' 'tcheb' 'Static polynomial nonlinearity representation'} ...
         {'threshold' 10^-10 'Threshold on SSE for terminitiaon of the iterative search'} ...
         {'decimate_option' 'iir' 'Available options are: no, ideal, iir'} ...
         {'use_vel_zeroth_exp' 'yes' 'Use the zero-th order basis expansion of delayed velocity: no, yes'} ...
         {'static_nl_init' [] 'Initial values for static NL parameters: empty (to be set by the algorithm), vector of last estimates'} ...
         {'static_nl_param_init' 'unityGain' 'Static NL initialization. Available options are: random, hwrFitKJ, hwrFitEST, unityGain'} ... % On Sept. 4, 2015, I also added the initialization option of 'unityGain'
     };
 if arg_parse(options,varargin);
     return
 end
% Author: Ehsan Sobhani  
% Date: August, 29th, 2015 Ver 0.0

%% Modifications:
%-- Dec. 05, 2015: I added a Least Squares Solution with the Regressor (i.e., identifying the entire parameters teta = alpha*c, without estimating the individual elements alpha and c).
%-- Dec. 17, 2015: I augmented the output structures for LPV static NL and LPV Laguerre-IRF with additional fields required for simulating the LPV Hammerstein system's response 
%-- May  05, 2021: I applied changes to the naming of properties and function arguments to make them consistent with the latest version of the NLID toolbox
%                  For example, nldatObj.data to nldatObj.dataSet
%% Initialization
it = max_iter; %1000;
ts = get(z,'domainIncr');
fs = 1/ts;
decimated_ts = ts*decimation_ratio;
decimated_fs = fs/decimation_ratio;

%-- Demeaning the input/output data
z = z-mean(z);
input = z(:,1);
output = z(:,2);

%% Data Decimation for Identification
input = decimate_kian(input,decimation_ratio,decimate_option);
output = decimate_kian(output,decimation_ratio,decimate_option);
rho = decimate_kian(rho,decimation_ratio,decimate_option);
nsamp = size(rho.dataSet,1);

%% Normalizing the input of the static NL of the Hammerstein system 
avg_i = (max(input.dataSet)+min(input.dataSet))/2;
rng_i = max(input.dataSet) - min(input.dataSet);
% avg_i = 0;
% rng_i = 5;
input_n = (input.dataSet - avg_i)*2/rng_i;

%% Normalizing scheduling variable
avg_rho = (max(rho.dataSet)+min(rho.dataSet))/2;
rng_rho = max(rho.dataSet) - min(rho.dataSet);
% avg_rho = -0.125;
% rng_rho = +0.75;
rho_n = (rho.dataSet - avg_rho)*2/rng_rho;

%% Computing the lagged, normalized SV and expanded input for the LPV Hammerstein system
Lr = round(irf_len_r / (decimated_ts));
if nside_r == 1
    lags_r = decimated_ts * (0:Lr);
else
    lags_r = decimated_ts * (-Lr:1:Lr);
end

nLags_r = length(lags_r);

rho_n_lags = zeros(nsamp,length(nLags_r));
input_n_lags = zeros(nsamp,length(nLags_r));
for d = 1:nLags_r
    rho_n_lags(:,d) = del(rho_n,decimated_fs,lags_r(d));
    input_n_lags(:,d) = del(input_n,decimated_fs,lags_r(d));
end

%% Computing the Laguerre Basis Functions up to order q
f = laguerre(nLags_r,q,alfa);

%% Initializing and Computing the Tchebychev expansions of:
%++ A) The normalized delayed velocity and its lags to order n
% Gn = zeros(nsamp,n+1,nLags_r);
% for d = 1:nLags_r
%     Gn(:,:,d) = multi_tcheb(input_n_lags(:,d),n);
% end

switch use_vel_zeroth_exp
    case 'no'
        Gn = zeros(nsamp,n,nLags_r);
    case 'yes'
        Gn = zeros(nsamp,n+1,nLags_r);
end

for d = 1:nLags_r
    temp = multi_tcheb(input_n_lags(:,d),n);
    switch use_vel_zeroth_exp
        case 'no'
            Gn(:,:,d) = temp(:,2:end);                 %-- we do not use the constant (zero-th order) term for delayed velocity expansion
        case 'yes'
            Gn(:,:,d) = temp;
    end
end

%++ B) The normalized Scheduling Variable and its lags to order p 
Gp = zeros(nsamp,p+1,nLags_r);
for d = 1:nLags_r
    Gp(:,:,d) = multi_tcheb(rho_n_lags(:,d),p);
end

%++ C) The normalized Scheduling Variable to order m 
Gm =  multi_tcheb(rho_n,m);

%% Computing the Z vector of expanding the reflex pathway (see my notes for formulas)
%-- On Sept. 22nd, I replaced all n with nF and conditioned nF on use_vel_zeroth_exp
switch use_vel_zeroth_exp
    case 'no'
        nF = n;
    case 'yes'
        nF = n+1;
end

Z = zeros(nsamp,(q+1)*(nF)*(p+1));

Gn_delay_pages = zeros(nsamp,nLags_r);
Gp_delay_pages = zeros(nsamp,nLags_r);

cnt = 0;
for cnt_n = 1:nF      
    for cnt_p = 1:p+1
        for cnt_q = 1:q+1
            cnt = cnt+1;
            for d = 1:nLags_r
                Gn_delay_pages(:,d) = Gn(:,cnt_n,d);
                Gp_delay_pages(:,d) = Gp(:,cnt_p,d);
            end
            GnGp = Gn_delay_pages.*Gp_delay_pages;
            Z(:,cnt) = GnGp * f(:,cnt_q);
        end
    end
end

%% The Regressor of the LPV-Hammerstein
Phi_lpvHamm_expand = zeros(nsamp,(q+1)*(nF)*(p+1)*(m+1));

for k = 1:nsamp
    Phi_lpvHamm_expand(k,:) = kron(Z(k,:),Gm(k,:));
end
    
% col = 0;
% for d = 1:nLags_r
%     for nc = 1:n+1   %-- nc is the counter for n
%         for mc = 1:m+1
%             for prc = 1:p_r+1
%                 col = col + 1;
%                 Phi_lpvHamm_expand(:,col) = Gn(:,nc,d).*Gm(:,mc,d).*Gpr(:,prc);
%             end
%         end
%     end
% end

%% Setting the output variable
o = output.dataSet;
tqI_r = o;

%% 
% disp('LS with All Regressor ...')
% teta_hat = lscov(Phi_lpvHamm_expand,tqI_r);
% tqI_r_hat_total = Phi_lpvHamm_expand * teta_hat;
% vafTotal = vaf(tqI_r,tqI_r_hat_total)

%% Iterative Separable Least Squares for identification of the reflex pathway
%% Decomposition of the parameters of the LPV static nonlinearity, alpha [n-by-(p+1)] 
%% and the parameters of the LPV Laguerre IRF, c [(q+1)-by-(m+1)]
alpha_hat = zeros(nF*(p+1),it);
alpha_hat_std = zeros(nF*(p+1),it);
 
c_hat = zeros((q+1)*(m+1),it);
c_hat_std = zeros((q+1)*(m+1),it);

if isempty(static_nl_init)
    switch static_nl_param_init
        case 'random'
            alpha0 = randn(nF*(p+1),1);
        case 'hwrFitKJ'
            %++ Kian's function for polynomial (chebychev) fit to HWR 
            alphaHWR = halfwave_rectifier_tchebychev(min(input.dataSet),max(input.dataSet),nF-1);   %--- This function generates (nF-1)+1 parameters!
            alpha0 = repmat(alphaHWR,(p+1),1);
        case 'hwrFitEST'
            %++ Ehsan's function (written on Spet. 4th, 2015) for polynomial (chebychev) fit to HWR 
            % slope = 1; threshold = 0;
            slope = 1; threshold = 0;
            plotFlag = 0; %1;
            alphaHWR = hwr_fit(min(input.dataSet),max(input.dataSet),slope,threshold,nF-1,'chebychev',plotFlag);
            % alphaHWR = alphaHWR(2:end);  %-- Excluding the constant term
            alpha0 = repmat(alphaHWR,(p+1),1);
        case 'unityGain'  %-- See my notes for reasons and details
            % alpha0 = zeros(nF*(p+1),1);
            % alpha0(1) = 1;
            alpha0 = zeros(nF,1); alpha0(2) = 1;
            alpha0 = repmat(alpha0,(p+1),1);
        otherwise 
            disp('This case is not supported!')
            return;
    end
else
    alpha0 = static_nl_init; 
end
%-- Initialization of the iterative algorithm
alpha0 = alpha0/norm(alpha0);
s1 = 10^10;
s2 = 10^10;
%-- Start of the iterative algorithm
for i = 1:it
    if i==1
        alpha0 = alpha0 .* wNL;  %-- Making the redundant alpha parameters zero in the construction of the Gamma_alpha matrix
        Gamma_alpha = Gamma_alpha_construct_laguerre(alpha0,q+1,m+1,nF,p+1);
        Phi_alpha_it = Phi_lpvHamm_expand*Gamma_alpha;
        %++ Using only the non-zero weighted elements of the regressor 
        Phi_alpha_it_w = zeros(nsamp,sum(wIRF));
        cntIRF = 0;
        for ii = 1:length(wIRF)
            if wIRF(ii) == 1
                cntIRF = cntIRF + 1;
                Phi_alpha_it_w(:,cntIRF) =  Phi_alpha_it(:,ii);
            end
        end
        [chat,chatSTD] = lscov(Phi_alpha_it_w,tqI_r);
        c_hat(wIRF==1,i) = chat;
        c_hat_std(wIRF==1,i) = chatSTD;
        % sse_c = tqI_r'*tqI_r - c_hat(:,i)'*Phi_alpha_it_w'*tqI_r;
        sse_c = tqI_r'*tqI_r - chat'*Phi_alpha_it_w'*tqI_r;
        
    else
        alpha_hat(:,i-1) = alpha_hat(:,i-1) .* wNL;  %-- Making the redundant alpha parameters zero in the construction of the Gamma_alpha matrix
        Gamma_alpha = Gamma_alpha_construct_laguerre(alpha_hat(:,i-1),q+1,m+1,nF,p+1);
        Phi_alpha_it = Phi_lpvHamm_expand*Gamma_alpha;
        %++ Using only the non-zero weighted elements of the regressor 
        Phi_alpha_it_w = zeros(nsamp,sum(wIRF));
        cntIRF = 0;
        for ii = 1:length(wIRF)
            if wIRF(ii) == 1
                cntIRF = cntIRF + 1;
                Phi_alpha_it_w(:,cntIRF) =  Phi_alpha_it(:,ii);
            end
        end
        
        [chat,chatSTD] = lscov(Phi_alpha_it_w,tqI_r);
        c_hat(wIRF==1,i) = chat;
        c_hat_std(wIRF==1,i) = chatSTD;
        % sse_c = tqI_r'*tqI_r - c_hat(:,i)'*Phi_alpha_it_w'*tqI_r;
        sse_c = tqI_r'*tqI_r - chat'*Phi_alpha_it_w'*tqI_r;
        
    end
    
    c_hat(:,i) = c_hat(:,i) .* wIRF;  %-- Making the redundant alpha parameters zero in the construction of the Gamma_c matrix
    Gamma_c = Gamma_c_construct_laguerre(c_hat(:,i),q+1,m+1,nF,p+1);
    Phi_c_it = Phi_lpvHamm_expand*Gamma_c;
    
    %++ Using only the non-zero weighted elements of the regressor 
    Phi_c_it_w = zeros(nsamp,sum(wNL));
    cntNL = 0;
    for ii = 1:length(wNL)
        if wNL(ii) == 1
            cntNL = cntNL + 1;
            Phi_c_it_w(:,cntNL) =  Phi_c_it(:,ii);
        end
    end
    
    [ahat,ahatSTD] = lscov(Phi_c_it_w,tqI_r);
    alpha_hat(wNL==1,i) = ahat;
    alpha_hat_std(wNL==1,i) = ahatSTD;
    % sse_alpha = tqI_r'*tqI_r - alpha_hat(:,i)'*Phi_c_it_w'*tqI_r;
    sse_alpha = tqI_r'*tqI_r - ahat'*Phi_c_it_w'*tqI_r;
    
    %==== Normalization step
    %-- Finding non-zero indices
    idxNL = find(wNL==1);
    idxIRF = find(wIRF==1);
    
    % h = sign(alpha_hat(1,i));        %-- This always selects the 1st one
    h = sign(alpha_hat(idxNL(1),i));   %-- This selects the 1st non-zero element
    
    % c_hat(:,i) = h * c_hat(:,i)*norm(alpha_hat(1:nF*(p+1),i));
    % alpha_hat(:,i) = alpha_hat(:,i) / norm(alpha_hat(1:nF*(p+1),i)) * h;
    c_hat(idxIRF,i) = h * c_hat(idxIRF,i)*norm(alpha_hat(idxNL,i));
    alpha_hat(idxNL,i) = alpha_hat(idxNL,i) / norm(alpha_hat(idxNL,i)) * h;
    if (s1-sse_c < threshold) && (s2-sse_alpha < threshold)
        break
    end
    s1 = sse_c;
    s2 = sse_alpha;
end

it = i;
% disp(['Terminated at iteration ',num2str(it)]);
fIter = it;

alpha = alpha_hat(:,it);
alpha_std = alpha_hat_std(:,it);
c = c_hat(:,it);
c_std = c_hat_std(:,it);

alpha_matrix = zeros(nF,p+1);
for cnt_n = 1:nF
    ind1 = (cnt_n-1)*(p+1) + 1;
    ind2 = cnt_n*(p+1);
    alpha_matrix(cnt_n,:) = alpha(ind1:ind2,1)';
end

c_matrix = zeros(m+1,q+1);
for cnt_q = 1:q+1 
    ind1 = (cnt_q-1)*(m+1) + 1; 
    ind2 = cnt_q*(m+1); 
    c_matrix(:,cnt_q) = c(ind1:ind2,1); 
end

static_nl.svNormalization = [avg_rho,rng_rho];
static_nl.data = alpha;
static_nl.dataSTD = alpha_std;
static_nl.coeffs = alpha_matrix;
static_nl.useZerothInpExp = use_vel_zeroth_exp;
static_nl.inputExpOrder = n;
static_nl.inputNormalization = [avg_i,rng_i];
static_nl.svExpOrder = p;

h_r.svNormalization = [avg_rho,rng_rho];
h_r.data = c;
h_r.dataSTD = c_std;
h_r.coeffs = c_matrix;
h_r.LaguerreAlfa = alfa;
h_r.LaguerreExpOrder = q;
h_r.svExpOrder = m;
h_r.nLags = nLags_r;
h_r.Ts = decimated_ts;

yp = Phi_c_it * alpha;
yp = nldat(yp,'domainIncr',decimated_ts);


