function y = pvnlSim(static_nl,u,rho)
%% Ehsan Sobhnai, May 12, 2021, extracted from the function pvHammLaguerreSim 
%++ This function simulates the output of a static PVNL model in response to the input u and scheduling variable rho

%-- Inputs: 
%       (1) static_nl: The structure containing the PVNL model
%       (2) u: input signal (e.g., perturbation velocity)
%       (3) rho: scheduling variable (e.g., joint position) 

%-- Output:
%       (1) y: predicted response of the model to input u and SV of rho

%% Modifications:
%-- May  05, 2021: I applied changes to the naming of properties and function arguments to make them consistent with the latest version of the NLID toolbox
%                  For example, nldatObj.data to nldatObj.dataSet

%% Body of the code
%-- Reading the structural blocks of the LPV Hammerstein system (LPV NL + LPV IRF)
nSamples = length(u.dataSet);
ts = u.domainIncr;

%% Reading the static PVNL parameters and simulating it output, z
avg_rho = static_nl.svNormalization(1);
rng_rho = static_nl.svNormalization(2);
avg_i = static_nl.inputNormalization(1);
rng_i = static_nl.inputNormalization(2);
alpha = static_nl.coeffs;
n = static_nl.inputExpOrder;
p = static_nl.svExpOrder;

%-- The parameters of static PVNL in matrix alpha must be ordered as [inputExpOrder,svExpOrder]
[inputExpOrder,svExpOrder] = size(alpha);
if svExpOrder ~=p+1
    disp('The expansion order of SV does not match the number of columns in alpha!')
    y = [];
end

switch static_nl.useZerothInpExp
    case 'no'
        if inputExpOrder ~= n
            disp('The expansion order of input does not match the number of rows in alpha!')
            y = [];
            return;
        end
    case 'yes'
        if inputExpOrder ~= n+1
            disp('The expansion order of input does not match the number of rows in alpha!')
            y = [];
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

y = nldat(z,'domainIncr',ts);

