function [yp,zp] = pvHammLaguerreSim(model,u,rho)
%% Ehsan Sobhnai, May 12, 2021, just renamed the function created Feb. 26, 2021 from lpvHammLaguerreSim to pvHammLaguerreSim 
%++ This function simulates the output of a LPV Hammerstein Laguerre model in response to the input u and scheduling variable rho

%-- Inputs: 
%       (1) model: The object/structure containing the model
%       (2) u: input signal (e.g., perturbation velocity
%       (3) rho: scheduling variable 

%-- Output:
%       (1) yp: predicted response of the model to input u and SV of rho

%% Modifications:
%-- May 05, 2021: I applied changes to the naming of properties and function arguments to make them consistent with the latest version of the NLID toolbox
%                  For example, nldatObj.data to nldatObj.dataSet
%-- May 12, 2021: I added a new output to the model which is the predicted output of the static NL block, named zp 

%% Body of the code
%-- Reading the structural blocks of the LPV Hammerstein system (LPV NL + LPV IRF)
static_nl = model.static_nl;
h_r = model.h_r;

nSamples = length(u.dataSet);

%+++ Few sanity checks on the model
if static_nl.svNormalization ~= h_r.svNormalization
    disp('The SV Normalization Factors for Static NL and IRF-Laguerre are not the same!')
    yp = [];
    return;
else
    avg_rho = static_nl.svNormalization(1);
    rng_rho = static_nl.svNormalization(2);
end

ts = h_r.Ts;
fs = 1/ts;
if ts ~= u.domainIncr
    disp('The sampling time of the input does not match that of the IRF dynamics of the Hammerstein cascade!')
    yp = [];
    return;
end

%% Reading the LPV Static NL parameters and simulating it output, z
avg_i = static_nl.inputNormalization(1);
rng_i = static_nl.inputNormalization(2);
alpha = static_nl.coeffs;
n = static_nl.inputExpOrder;
p = static_nl.svExpOrder;

%-- The parameters of LPV static NL in matrix alpha must be ordered as [inputExpOrder,svExpOrder]
[inputExpOrder,svExpOrder] = size(alpha);
if svExpOrder ~=p+1
    disp('The expansion order of SV does not match the number of columns in alpha!')
    yp = [];
end

switch static_nl.useZerothInpExp
    case 'no'
        if inputExpOrder ~= n
            disp('The expansion order of input does not match the number of rows in alpha!')
            yp = [];
            return;
        end
    case 'yes'
        if inputExpOrder ~= n+1
            disp('The expansion order of input does not match the number of rows in alpha!')
            yp = [];
            return;
        end
end

%++ Normalizing the input and scheduling variable
u_n = (u.dataSet - avg_i)*2/rng_i;
rho_n = (rho.dataSet - avg_rho)*2/rng_rho;

%++ Chebychev Basis Exapnsions of normalized Input and SV
Gn = multi_tcheb(u_n,n); 
if strcmp(static_nl.useZerothInpExp,'no')
    Gn = Gn(:,2:end);
end

Gp = multi_tcheb(rho_n,p);

%++ Generating the output of the LPV static NL, z
z = zeros(nSamples,1);
for k = 1:nSamples
    sumZ = 0;
    for i = 1:size(Gn,2)
        for j = 1:size(Gp,2)
            sumZ = sumZ + alpha(i,j)*Gn(k,i)*Gp(k,j); 
        end
    end
    z(k,1) = sumZ;
end

zp = nldat(z,'domainIncr',ts);
zp.comment = 'Predicted output of the static PVNL of the Hammerstein cascade'; 

%% Reading the LPV Laguerre Dynamics parameters and simulating it output, i.e. system output yp, using z as its input.
c = h_r.coeffs;
alfa = h_r.LaguerreAlfa;
q = h_r.LaguerreExpOrder;
m = h_r.svExpOrder;
nLags = h_r.nLags;
  
%-- Generating the Laguerre Basis Functions up to order q
f = laguerre(nLags,q,alfa);

%-- Generating Chebychev Basis Expansions of SV
Gm = multi_tcheb(rho_n,m);

%-- Computing the lagged inputs of the LPV Hammerstein dynamics, which are the lagged outputs of the LPV static NL.  
lags = h_r.Ts * (0:nLags);
z_lags = zeros(nSamples,length(nLags));

for d = 1:nLags
    z_lags(:,d) = del(z,fs,lags(d));
end

%-- Computing the Zl[nsamp,(q+1)] vector of covolution sum of Laguerre bases and lagged outputs of the LPV static NL. 
ZL = zeros(nSamples,q+1);  %-- I named it with Capital ZL (output of static NL convolved with Laguerre)

for k = 1:nSamples
    for cnt_q = 1:q+1
        ZL(k,cnt_q) = f(:,cnt_q)'*z_lags(k,:)';
    end
end

%% Computing the output of the LPV Hammerstein system
yp = zeros(nSamples,1);
for k = 1:nSamples
    sumY = 0;
    for cnt_m = 1:m+1
        for cnt_q = 1:q+1
            sumY = sumY + c(cnt_m,cnt_q)*Gm(k,cnt_m)*ZL(k,cnt_q);
        end
    end
    yp(k,1) = sumY;
end

% yp = h_r.Ts * yp;

yp = nldat(yp,'domainIncr',ts);
yp.comment = 'Predicted output of the PV Hammerstein cascade'; 

