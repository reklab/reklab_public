classdef pvpc < pvm
    % pvpc - parameter-varying parallel-cascade model of joint stiffness for NLID toolbox.
    
    % pvpc model has two pathways: 
    % Intrinsic pathway represnted by a PV IRF model (pvirf class)
    % Reflex pathway represented by a PV Hammerstein model (pvnlbl class) 
    % These will be stored in the elements property as {pvirf;pvnlbl}.
    
    % Identification methods supported:
    % 1) NPPVPC: Non-Parametric Parameter-Varying Parallel-Cascade method
    %    See E. Sobhani Tehrani and R.E. Kearney, "Non-Parametric Nonlinear Parameter-Varying Parallel-Cascade Identification of Dynamic Joint Stiffness"
    %        submitted to IEEE TBME, 11 Feb. 2022
    
    properties
        identData = NaN;     %-- I/SV/O data used for identification must be kept within object. 
        identTQt = NaN;      %-- remove this as related to data set used for ID. --> The identified output from model identification/training must be kept within object. 
        identVAF = NaN;      %-- remove this as related to data set used for ID. --> The performance indicator of model identification/training must be kept within object.
        identTQi = NaN;
        identTQr = NaN;
    end
    
    methods
        function sys = pvpc(a,varargin)
            sys.notes = 'Non-Parametric Parameter-Varying Parallel-Cascade (PV-PC) model of joint stiffness';
            sys.comment = 'PC stiffness model elements are modulate by a scheduling variable (SV)';
            
            sys.parameterSet(1)=param('paramName','idMethod', ...
                'paramDefault','',...
                'paramHelp', 'A method that uses orthogonal projections and an iterative method to identify PV models of intrinsic and reflex stiffness',...
                'paramType','select',...
                'paramLimits', {'','nppv-pc','sspv-pc'});
            
            sys.elements = {pvirf; pvnlbl};
            
            if nargin==0
                return
            elseif isa(a,'pvpc')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end
        end
        
        
        %% Setting system parameters specific to each identification method
        %++ Set method first so the parameters are properly defined
        function sys = set(sys, varargin)
            %-- First, extract the name of the system using MATLAB's inherent function inputname
            sysname = inputname(1); 
            %-- Find the identification method from the parameters cell array
            iMethod = find(strcmp(varargin,'idMethod'));
            if ~isempty(iMethod)
                if length(varargin)>=iMethod+1 %-- I am not quite sure what this If condition is meant for
                    %-- Call the setidMethod to set the parameters specific
                    %-- to the selected method.
                    sys = setidMethod(sys,varargin{iMethod+1});
                else
                    sys.parameterSet = setval(sys.parameterSet,varargin);
                end
            end
            set@nltop(sys, varargin{:});
            assignin('caller',sysname,sys);
        end
        
        %% Function to set parameters specific to each identification method
        function sys = setidMethod(sys,newMethod)          
            methodList= {'','nppv-pc','sspv-pc'};
            %++ ToDo: The above line must be replaced reading the limits from sys object
            %++       something like this:
            % P = sys.parameterSet;
            % methodList = P(1).paramLimits; %-- In principle, we have to read the in index of idMethod instead of assuming it is 1
            if ~any(strcmp(newMethod, methodList)) % isempty(newMethod) || ~any(strcmp(newMethod, methodList))
                disp('Not a valid identification method option');
                disp (['Available identification method options are: [' strjoin(methodList) ']' ]);
                return
            end
           
            switch newMethod
                case ''
                    %--- no need to add any parameters as ident is not defined
                case 'nppv-pc'
                    j = length(sys.parameterSet);
                    %++ Parameters to identify intrinsic pathway
                    sys.parameterSet(j+1) = param('paramName','irf_len_i',...
                        'paramDefault',0.04,... 
                        'paramHelp','Length of intrinsic IRF in seconds',...
                        'paramType','number',...
                        'paramLimits',[0 0.1]);
                    
                    sys.parameterSet(j+2) = param('paramName','nside_i',...
                        'paramDefault',2,...
                        'paramHelp','Number of sides of the intrinsic IRF',...
                        'paramType','number',...
                        'paramLimits',{1, 2});
                    
                    sys.parameterSet(j+3) = param('paramName','p_i',...
                        'paramDefault',2,... 
                        'paramHelp','Order of the intrinsic IRF expansion w.r.t. SV',...
                        'paramType','number',...
                        'paramLimits',[0 11]);
                    
                    %++ Parameters to identify intrinsic pathway
                    sys.parameterSet(j+4) = param('paramName','irf_len_r',...
                        'paramDefault',0.8,... %0.6,...
                        'paramHelp','Length of reflex IRF in seconds',...
                        'paramType','number',...
                        'paramLimits',[0 inf]);
                    
                    sys.parameterSet(j+5) = param('paramName','nside_r',...
                        'paramDefault',1,...
                        'paramHelp','Number of sides of the reflex IRF',...
                        'paramType','number',...
                        'paramLimits',{1, 2});
                                      
                    %-- Change name from alfa to laguerre_decay
                    sys.parameterSet(j+6) = param('paramName','alfa',...
                        'paramDefault',0.54,... %0.8,...
                        'paramHelp','Constant parameter of the Laguerre expansion',...
                        'paramType','number',...
                        'paramLimits',[0.01 1]); %-- The theoretical lower limit is any number bigger than 0
                    
                    %-- The theoretical feasible range for Laguerre exapansion order is [1, inf].
                    %-- But any expansion order < 5 is pratically useless.
                    %-- Any expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from q to laguerre_order
                    sys.parameterSet(j+7) = param('paramName','q',...
                        'paramDefault',8,... %10,...
                        'paramHelp','Order of the Laguerre expansion of the IRF',...
                        'paramType','number',...
                        'paramLimits',[5 16]); 
                                     
                    %-- The theoretical feasible range for SV exapansion order of Laguerre coefficients is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from m to laguerre_sv_order
                    sys.parameterSet(j+8) = param('paramName','m',...
                        'paramDefault',8,... %3,...
                        'paramHelp','Order of SV expansion for IRF Laguerre coefficients',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 corresponds to time-invariant IRF
                    
                    %-- The theoretical feasible range for Input exapansion order of static nl is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from n to nl_input_order
                    sys.parameterSet(j+9) = param('paramName','n',...
                        'paramDefault',9,... %6,...
                        'paramHelp','Order of input expansion for static nonlinearity',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 converts Hammerstein to linear dynamics (TI/PV IRF)
                    
                    %-- The theoretical feasible range for SV exapansion order of static NL coefficients is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from p to nl_sv_order
                    sys.parameterSet(j+10) = param('paramName','p',...
                        'paramDefault',7,... %2,...
                        'paramHelp','Order of SV expansion for static nonlinearity coefficients',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 corresponds to time-invariant IRF
                    
                    sys.parameterSet(j+11) = param('paramName','max_iter',...
                        'paramDefault',500,...
                        'paramHelp','Maximum number of iterations',...
                        'paramType','number',...
                        'paramLimits',[1 1000]);
                    
                    sys.parameterSet(j+12) = param('paramName','threshold',...
                        'paramDefault',10^-10,...
                        'paramHelp','Threshold on SSE for terminitiaon of the iterative search',...
                        'paramType','number',...
                        'paramLimits',[0 1]);
                    
                case 'sspv-pc'
                    sys.parameterSet = param;
                    j = length(sys.parameterSet);
                    sys.parameterSet(j+1)=param('paramName','idMethod', ...
                        'paramDefault','sspv-pc',...
                        'paramHelp', 'This is subspace method to identify the PV-PC model of stiffness with TI reflex dynamics');
                    disp('The code for the SSPV-PC method has not been integrated into NLID toolbox yet.');
                    
                otherwise
                    %-- If other ident methods are added, put their specific parameters as a new switch case here.    
            end
        end
        
        %% Implementation of various PV-PC identification methods 
        % For example, based on a comment from meeting with Rob, we need something like: nlident(sys,z,'method','npv_hamm')
        function sys = nlident(sys, z, varargin)
            sys.identData = z;
            
            iMethod = find(strcmp(varargin,'idMethod'));
            newMethod = varargin(iMethod+1);
            %++ ToDo: The above line must be replaced reading the limits from sys object
            %++ something like this:
            % P = sys.parameterSet;
            % methodList = P(1).paramLimits; %-- In principle, we have to read the in index of idMethod instead of assuming it is 1
            methodList= {'nppv-pc','sspv-pc'};
            if isempty(newMethod) || ~any(strcmp(newMethod, methodList))
                disp('Not a valid identification method option');
                disp (['Available identification method options are: [' strjoin(methodList) ']' ]);
                return
            end
            
            %== Decimation should naturaly be applied to input triplet data (input, sv, output) - TBD 
            iDecimation = find(strcmp(varargin,'decimation'));
            decimation = varargin(iDecimation+1); decimation = decimation{1,1};
            
            %== Reflex delay (rDelay) should naturaly become an attribute of the model or a (hyper)parameter of the identification method - TBD  
            iDelay = find(strcmp(varargin,'rDelay'));
            rDelay = varargin(iDelay+1); rDelay = rDelay{1,1};
                       
            switch newMethod{1,1}             
                case 'nppv-pc'
                    io = z(:,[1,3]); %-- Input is the 1st and Output is the last data
                    rho = z(:,2);    %-- schedVar is the 2nd data
                    %++ ToDo: The above lines must be removed by updating
                    %         the identification function to accept 
                    %         triplet [input,sv,output] instead of rho as a
                    %         separate signal.
                    
                    %++ ToDo: The following lines shall not be needed 
                    %         when "hammerstein_lpv_laguerre_v02" is
                    %         updated to accept the sys object.
                    
                    %++ For now, extract identification parameters from the sys object outside the function
                    %++ 1st: Extract INTRINSIC stiffness paramaters
                    irf_len_i = get(sys,'irf_len_i');  %-- in seconds
                    nside_i = get(sys,'nside_i');
                    p_i = get(sys,'p_i');
                    
                    %++ 2nd: Extract REFLEX stiffness paramaters
                    irf_len_r = get(sys,'irf_len_r');
                    alfa = get(sys,'alfa');
                    q = get(sys,'q');
                    m = get(sys,'m');
                    n = get(sys,'n');
                    p = get(sys,'p'); 
                    max_iter = get(sys,'max_iter');   
                    
                    %++ ToDo: We have to decide which of the following 
                    %         parameters must be added to object or not
                    %== The weight functions of PVIRF and PVNL to allow
                    %   making specific parameters of IRF or NL independent from SV
                    wIRF = ones((m+1)*(q+1),1);
                    wNL = ones((n+1)*p,1);
                    
                    %== This we need to add to ident parameters as an initialization value for the NACS iteration of NPNPV-H 
                    static_nl_param_init = 'hwrFitEST';
                    %== The following two parameters are actually not used.
                    %   Delete them later from "hammerstein_lpv_laguerre_v02" inputs
                    use_vel_zeroth_exp = 'yes';
                    % static_nl_init = []; 
                    
                    [h_i,h_r,static_nl,TQIhat,TQRhat,TQRProj] = pcasc_lpv_laguerre(io,rho,wNL,wIRF,'irf_len_r',irf_len_r,'irf_len_i',irf_len_i,'n',n,'p',p-1,'p_i',p_i,'m',m,'q',q,'alfa',alfa,'delay',rDelay,...
                                                                             'static_nl_param_init',static_nl_param_init,'use_vel_zeroth_exp',use_vel_zeroth_exp,'nside_i',nside_i,'max_iter',max_iter);
                    
                    %% Set output predictions and VAF(%) for the identified model and 
                    sys.identTQi = TQIhat;
                    sys.identTQr = TQRhat;
                    sys.identTQt = TQIhat + TQRhat;
                    %++ Calculate VAF
                    o_d = decimate_kian(z(:,3),decimation);
                    v = vaf(o_d,sys.identTQt);
                    sys.identVAF = v.dataSet;
                    
                    %% Setup the objects based on the identification results
                    rho_d = decimate_kian(rho,decimation);
                    min_sv = min(rho_d.dataSet); max_sv = max(rho_d.dataSet);
                    
                    %++ INTRINSIC Pathway:
                    PVIRF_i_BASIS = mimobasis;
                    nside_i = get(sys,'nside_i');
                    if nside_i == 1
                        min_u = 0;
                        max_u = h_i.nLags - 1;
                    elseif nside_i == 2
                        min_u = - fix(h_i.nLags / nside_i);
                        max_u = + fix(h_i.nLags / nside_i);
                    else
                        disp('The number of sides must be either 1 or 2! Please correct and rerun!')
                        return;
                    end
                    PVIRF_i_BASIS = set(PVIRF_i_BASIS,'coeffs',h_i.data,...
                                      'svRange',[min_sv,max_sv],...
                                      'inputRange',[min_u,max_u],...
                                      'svExpType','tcheb',...
                                      'inputExpType','none',...
                                      'svExpOrder',h_i.svExpOrder,...
                                      'inputExpOrder',NaN);
                                  
                    PVIRF_i = pvirf;
                    PVIRF_i.nSides = nside_i;
                    PVIRF_i.elements = PVIRF_i_BASIS;
                    
                    %++ REFLEX Pathway: Extracting SV Tchebychev polynomials from static_nl and setup PVNL 
                    sv_polynoms_pvnl = cell(static_nl.inputExpOrder+1,1);
                    
                    i_d = decimate_kian(z(:,1),decimation);
                    min_u = min(i_d.dataSet); max_u = max(i_d.dataSet);

                    for j = 1:length(sv_polynoms_pvnl)
                        sv_polynoms_pvnl{j,1} = polynom('polyCoef',static_nl.coeffs(j,:)',...
                                                   'polyOrder',static_nl.svExpOrder,...
                                                   'polyType','tcheb');                             
                    end
                            
                    PVNL_BASIS = mimobasis;
                    PVNL_BASIS = set(PVNL_BASIS,'coeffs',sv_polynoms_pvnl,...
                                'svRange',[min_sv,max_sv],...
                                'inputRange',[min_u,max_u],...
                                'svExpType','tcheb',...
                                'inputExpType','tcheb',...
                                'svExpOrder',static_nl.svExpOrder,...
                                'inputExpOrder',static_nl.inputExpOrder);
                    
                    PVNL = pvnl;
                    PVNL.elements = PVNL_BASIS;
                    
                    %++ Extracting SV Tchebychev polynomials from h_r and setup PVIRF 
                    sv_polynoms_pvirf = cell(h_r.LaguerreExpOrder+1,1);
                    for j = 1:length(sv_polynoms_pvirf)
                        sv_polynoms_pvirf{j,1} = polynom('polyCoef',h_r.coeffs(:,j),... 
                                                   'polyOrder',h_r.svExpOrder,...
                                                   'polyType','tcheb');                             
                    end
                    min_u = 0; max_u = (h_r.nLags - 1)*h_r.Ts;
                                       
                    PVIRF_BASIS = mimobasis;
                    PVIRF_BASIS = set(PVIRF_BASIS,'coeffs',sv_polynoms_pvirf,...
                                      'svRange',[min_sv,max_sv],...
                                      'inputRange',[min_u,max_u],...
                                      'svExpType','tcheb',...
                                      'inputExpType','laguerre',...
                                      'svExpOrder',h_r.svExpOrder,...
                                      'inputExpOrder',h_r.LaguerreExpOrder,...
                                      'alfa',h_r.LaguerreAlfa);
                    
                    PVIRF = pvirf;
                    PVIRF.elements = PVIRF_BASIS;
                    
                    PVH_r = pvnlbl;
                    set(PVH_r,'elements',{PVNL,PVIRF});
                    PVH_r.identData = TQRProj;
                    PVH_r.identOutput = sys.identTQr;
                    v = vaf(TQRProj,sys.identTQr);
                    PVH_r.identVAF = v.dataSet;
                    
                    %++ Set the estimated intrinsic model (PVIRF_i) and the estimated reflex model (PVH_r) as the elements of the identified PV-PC stiffness model
                    set(sys,'elements',{PVIRF_i; PVH_r});                    
                    
                case 'sspv-h'
                    disp('The code for SubSpace Hammerstein identification method is to be implemented from E. Sobhani Tehrani et al., IEEE EMBC, 2015.');
                
                otherwise
                    disp('This method does not exist for PV-Hammerstein systems.');
            end   

        end
        
        %% Function to plot the PV Hammerstein system (cascade of PV NL and PV IRF)
        function plot(sys,varargin)
            options={{'n_bins_input' 40 'number of bins for a grid on input'} ...
                     {'n_bins_sv' 40 'number of bins for a grid on SV'} ...
            };
            if arg_parse(options,varargin)
                return
            end
            
            %++ Extract elements of the cascade: PVNL and PVIRF
            PVNL = sys.elements{1,1};
            PVIRF = sys.elements{1,2};
            
            subplot(1,2,1)
            plot(PVNL,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv)
            subplot(1,2,2)
            plot(PVIRF,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv)
            
        end
        
        %% Simulation method for the PV Hammerstein systems
        function sys = nlsim(sys, z, varargin)
            %-- The simulation function
            u = z{1,1};
            rho = z{1,2};
            sys.yp = pvHammLaguerreSim(model,u,rho);
        end
        
    end %--> End of methods
end     %--> End of classdef

%% Required functions/subroutines for identification
%- 1. NPPV-PC Identification Method pcasc_lpv_laguerre 
%- 2. NPPV-H Identification Method (hammerstein_lpv_laguerre_v02)
%- 3. Laguerre polynomials or basis functions calculator
%- 4. Half wave rectifier fit for initialization of static nonlinearity
%- 5. Data matrix (Gamma_alpha) constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
%- 6. Data matrix (Gamma_c)     constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
%- 7. NPNPV-H simulator function

%% NPPV-PC identification method
function [h_i,h_r,static_nl,TQ_I,TQ_R,TQ_R_Proj,stopIter]= pcasc_lpv_laguerre(z,rho,wNL,wIRF,varargin)
% This function estimates an LPV Parallel-Cascade Model of Total Joint Stiffness.
% The model structure identified is an LPV Hammerstein (LPV Static NL + LPV Laguerre Expansion of IRF) for Reflex Pathway and an LPV IRF for Intrinsic Pathway.
% The inputs (in the NLDAT object z) are:
% -- 1) z: Input/Output data including: 
%          The Perturbation Poisition (in rad)
%          Joint Torque in Response to Perturbations (in Nm)
% -- 2) rho: The scheduling variable
% -- 3) wNL: Weights related to NL parameters: 1 if a parameter is required and 0 if a parameter must be excluded (added on Aug. 28, 2015 to accomodate NL parameter selection)
% -- 4) wIRF: Weights related to NL parameters: 1 if a parameter is required and 0 if a parameter must be excluded (added on Aug. 28, 2015 to accomodate Laguerre IRF parameter selection)

% The outputs are:
%--- 1) h_i: The structure for LPV IRF of the intrinsic pathway. The structure has the following parameters: 1)h_i.data; 2)h_i.domainIncr; 3)h_i.nSides. 
%       If 2-sided IRF is used, the dimension of h_i.data is [(2*Li+1)*p_i,1], where Li is the length of intrinsic IRF and p_i is the order of expansion of intrinsic IRF w.r.t. scheduling variable.

%--- 2) h_r: The structure for LPV Laguerre exapnsion of the reflex pathway IRF. The structure has the following parameters: 1)h_r.lpv_laguerre; 2)h_r.domainIncr
%       The dimension of h_r.data is [m*q,1], where q is the order of the Laguerre expansion of the reflex IRF and m is the expansion order of the reflex IRF Laguerre w.r.t. scheduling variable.

%--- 3) static_nl is the 2D static nonlinearity with dimension [n*p,1].
%       n is the order of expansion of delayed velocity and p is the order of expansion of the scheduling variable.

%--- 4) TQ_I is the estimated intrinsic torque 

%--- 5) TQ_R is the estimated reflex torque

%--- 6) TQ_R_Proj is the projected reflex torque (Output of the Orthogonal Projection step)

%--- 7) stopIter is the iteration at which the algorithm stopped (or converged)

options={{'decimation_ratio' 10 'ratio to decimate data bt'} ...
         {'irf_len_i' 0.04 'length of intrinsic IRF in s'} ...
         {'nside_i' 2 'number of sides of intrinsic IRF'} ...
         {'p_i' 3 'order of SV expansion for intrinsic IRF'} ...
         {'irf_len_r' 0.6 'length of reflex IRF in s'} ...
         {'nside_r' 1 'number of sides of intrinsic IRF'} ...
         {'q' 14 'order of the Laguerre expansion of the reflex IRF'} ...
         {'alfa' 0.8 'constant parameter of the Laguerre expansion'} ...
         {'m' 2 'order of SV expansion for reflex IRF Laguerre exapnsion'} ...
         {'n' 8 'order of velocity expansion for static nonlinearity'} ...
         {'p' 3 'order of SV expansion for static nonlinearity'} ...
         {'max_iter' 500 'maximum number of iterations for reflex pathway identification'} ...
         {'delay' 0.05 'Reflex delay'} ...
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
% Date: May, 8th, 2014 Ver 0.0

%% Modifications:
%++ Aug. 28, 2015: I changed the code to accomodate the use of wNL and wIRF in order to allow selection of parameters for Nonlinearity and IRF  
%                  The main changes were made in the initialization of alpha_hat, alpha_hat_std, c_hat, c_hat_std as well as the regression matrices Phi_alpha_it and Phi_c_it 

%++ Aug. 29, 2015: Also, since the lengths of the vectors alpha and c changes, I had to change and initialize them as cells  

%++ July 14, 2016: I updated this code based on hammerstein_lpv_laguerre_v02 code

%% Initialization
% it = max_iter; %1000;
ts = get(z,'domainIncr');
fs = 1/ts;
decimated_ts = ts*decimation_ratio;
decimated_fs = fs/decimation_ratio;

%-- Demeaning the input/output data
z = z-mean(z);
pos = z(:,1);
trq = z(:,2);

% Define velocity signal as the input of the reflex pathway.
vel = ddt(pos);
% acc = ddt(vel);

% Finding delayed velocity (because of the reflex delay)
dvel = del(vel.dataSet,fs,delay);
dvel = nldat(dvel,'domainIncr',ts);

%% Data Decimation for Identification
pos = decimate_kian(pos,decimation_ratio,decimate_option);
dvel = decimate_kian(dvel,decimation_ratio,decimate_option);
trq = decimate_kian(trq,decimation_ratio,decimate_option);
rho = decimate_kian(rho,decimation_ratio,decimate_option);
nsamp = size(rho.dataSet,1);

%% Computing the lagged positions for the Intrinsic IRF 
Li = round(irf_len_i / (decimated_ts));

if nside_i == 1
    lags_i = decimated_ts * (0:Li);
else
    lags_i = decimated_ts * (-Li:1:Li);
end

nLags_i = length(lags_i);
pos_lags = zeros(nsamp,length(nLags_i));
for d = 1:nLags_i
    pos_lags(:,d) = del(pos.dataSet,decimated_fs,lags_i(d));
end

%% Normalizing delayed velocity as input of the reflex pathway
avg_r = (max(dvel.dataSet)+min(dvel.dataSet))/2;
rng_r = max(dvel.dataSet) - min(dvel.dataSet);
dvel_n = (dvel.dataSet - avg_r)*2/rng_r;

%% Normalizing scheduling variable
avg_rho = (max(rho.dataSet)+min(rho.dataSet))/2;
rng_rho = max(rho.dataSet) - min(rho.dataSet);
rho_n = (rho.dataSet - avg_rho)*2/rng_rho;

%% Computing the Tchebychev expansions of scheduling for intrinsic pathway
Gpi =  multi_tcheb(rho_n,p_i);

%% The Regressor of the Intrinsic Pathway
Phi_intrinsic_expand = zeros(nsamp,(p_i+1)*(nLags_i));

col = 0;
for d = 1:nLags_i
    for pic = 1:p_i+1
        col = col + 1;
        Phi_intrinsic_expand(:,col) = pos_lags(:,d).*Gpi(:,pic);        
    end
end

%% Computing the lagged, normalized SV and expanded input for the reflex Pathway
Lr = round(irf_len_r / (decimated_ts));
if nside_r == 1
    lags_r = decimated_ts * (0:Lr);
else
    lags_r = decimated_ts * (-Lr:1:Lr);
end

nLags_r = length(lags_r);

rho_n_lags = zeros(nsamp,length(nLags_r));
dvel_n_lags = zeros(nsamp,length(nLags_r));
for d = 1:nLags_r
    rho_n_lags(:,d) = del(rho_n,decimated_fs,lags_r(d));
    dvel_n_lags(:,d) = del(dvel_n,decimated_fs,lags_r(d));
end

%% Computing the Laguerre Basis Functions up to order q
f = laguerre(nLags_r,q,alfa);

%% Initializing and Computing the Tchebychev expansions of:
%++ A) The normalized delayed velocity and its lags to order n
% Gn = zeros(nsamp,n+1,nLags_r);
% for d = 1:nLags_r
%     Gn(:,:,d) = multi_tcheb(dvel_n_lags(:,d),n);
% end

switch use_vel_zeroth_exp
    case 'no'
        Gn = zeros(nsamp,n,nLags_r);
    case 'yes'
        Gn = zeros(nsamp,n+1,nLags_r);
end

for d = 1:nLags_r
    temp = multi_tcheb(dvel_n_lags(:,d),n);
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

%% The Regressor of the Reflex Pathway (or the Regressor of the LPV-Hammerstein)
Phi_reflex_expand = zeros(nsamp,(q+1)*(nF)*(p+1)*(m+1));

for k = 1:nsamp
    Phi_reflex_expand(k,:) = kron(Z(k,:),Gm(k,:));
end
    
%% Setting the output variable
o = trq.dataSet;

%% Orthogonal Projection 
% tic;
intrinsic = intrinsic_estimator(Phi_intrinsic_expand,Phi_reflex_expand,o);
% toc;
tqI = Phi_intrinsic_expand*intrinsic;
tqI_r = o - tqI;
%-- Demeaning the reflex torque for identification of reflex pathway
tqI_r = tqI_r - mean(tqI_r); 

%% Assignment of function outputs
h_i.normalization = [avg_rho,rng_rho];
h_i.data = intrinsic;
h_i.nSides = nside_i;
h_i.nLags = nLags_i;
h_i.svExpOrder = p_i;
TQ_I = nldat(tqI,'domainIncr',decimated_ts,'comment','Estimated Intrinsic Torque');
TQ_R_Proj = nldat(tqI_r,'domainIncr',decimated_ts,'comment','Orthogonal Projection of Reflex Torque');

%% Iterative Normalized Convex Search for the identification of the NPLPV-H model of reflex pathway
%% Here is the part I replaced with a function call to hammerstein_lpv_laguerre_v02
%-- Constructing the input/output data for identification of reflex pathway
zR = cat(2,dvel,TQ_R_Proj);
decimationR = 1;             %-- The data is already decimated at this stage where we want to identify the reflex pathway. So, decimation_ration must be set to 1 because no more decimation is required
[static_nl,h_r,TQ_R,fIter] = hammerstein_lpv_laguerre_v02(zR,rho,wNL,wIRF,'irf_len_r',irf_len_r,'n',n,'p',p,'m',m,'q',q,'alfa',alfa,'max_iter',max_iter,...
                                                                 'static_nl_param_init',static_nl_param_init,'use_vel_zeroth_exp',use_vel_zeroth_exp,'decimation_ratio',decimationR,'static_nl_init',static_nl_init);
 
disp(['Terminated at iteration ',num2str(fIter)]);
stopIter = fIter;                                                             

end

%% 2. NPPV-H Identification Method
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

end


%% 2. Laguerre polynomials or basis functions calculator
%++ This function generates Laguerre orthonormal basis functions
%++ Author: Ehsan Sobhani (10 April 2014)
%++ This is based on Maremaleris book OR formula (11) of his paper titled:
%++ "Identification of Nonlinear Biological Systems Using Laguerre Expansions of Kernels", Annals of Biomed. Eng., vol. 21, pp. 573-589, 1993.
function b = laguerre(irf_len,max_order,alfa)
%++ The inputs are:
    % 1) irf_len (integer number of samples)
    % 2) max_order (integer). The maximum order of Laguerre expansion.
    % 3) alfa (0<alfa<1). Laguerre parameter. 

%++ The outputs are:
    % 1) b (real with size [irf_len,max_order+1])

%=== Initialization
b = zeros(irf_len,max_order+1);
    
for t = 1:irf_len
    for j = 0:max_order
        gain = alfa^((t-j)/2) * (1-alfa)^(0.5);
        summation = 0;
        for k = 0:j
            argument = (-1)^k * combination(t,k) * combination(j,k) * alfa^(j-k) * (1-alfa)^k;
            summation = summation + argument;
        end
        b(t,j+1) = gain * summation;
    end
end

end


%% 3. Half wave rectifier fit for initialization of static nonlinearity
function [alpha,yhat,u,Gu_n,y] = hwr_fit(uMin,uMax,slope,threshold,order,polynomType,plotFlag)

u = (uMin:0.0001:uMax)';
u_n = u/uMax;

y = slope*(max(u,threshold) - threshold);

switch polynomType
    case 'power'
        aPolynom = polyfit(u,y,order);
        yhat = polyval(aPolynom,u);
        alpha = aPolynom; %.coeff;
        Gu_n = '';
    case 'chebychev'
        Gu_n = multi_tcheb(u_n,order);
        alpha = lscov(Gu_n,y);
        yhat = Gu_n*alpha;
end

if plotFlag
    figure;
    plot(u,y,u,yhat,'r'); xlabel('Input'); ylabel('NL Output')
    legend('HWR','Fitted Polynom')
end

end

%% 4. Data matrix constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
function Gamma_alpha = Gamma_alpha_construct_laguerre(alpha,q,m,n,p)
Gamma_alpha = zeros(q*m*n*p,q*m);

cnt = 0;
ind2 = 0;
for i = 1:n
    for j = 1:p
        cnt = cnt+1;
        row_block = diag(alpha(cnt,1)*ones(q*m,1));
        ind1 = ind2+1;
        ind2 = cnt*q*m;
        Gamma_alpha(ind1:ind2,:) = row_block;
    end
end

end

%% 5. Data matrix constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
function Gamma_c = Gamma_c_construct_laguerre(c,q,m,n,p)
Gamma_c = zeros(q*m*n*p,n*p);

cnt = 0;
ind2 = 0;
for i = 1:n
    for j = 1:p
        cnt = cnt+1;
        ind1 = ind2 + 1;
        ind2 = cnt*q*m;
        Gamma_c(ind1:ind2,cnt) = c;
    end
end

end

%% 6. NPNPV-H simulator function
function yp = pvHammLaguerreSim(model,u,rho)
%% Ehsan Sobhnai, May 12, 2021, just renamed the function created Feb. 26, 2021 from lpvHammLaguerreSim to pvHammLaguerreSim 
%++ This function simulates the output of a LPV Hammerstein Laguerre model in response to the input u and scheduling variable rho

%-- Inputs: 
%       (1) model: The object/structure containing the model
%       (2) u: input signal (e.g., perturbation velocity
%       (3) rho: scheduling variable 

%-- Output:
%       (1) yp: predicted response of the model to input u and SV of rho

%% Modifications:
%-- May  05, 2021: I applied changes to the naming of properties and function arguments to make them consistent with the latest version of the NLID toolbox
%                  For example, nldatObj.data to nldatObj.dataSet

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

end