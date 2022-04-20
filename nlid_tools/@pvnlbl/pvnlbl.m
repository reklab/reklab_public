classdef pvnlbl < pvm
    % pvnlbl - parameter-varying nonlinear-linear model (Hammerstein model) for NLID toolbox.
    
    % pvnlbl model comprises a cascade of a pvnl model followed by a pvirf model 
    % stored in the elements property {pvnl,pvirf}.
    
    % Identification methods supported:
    % 1) Non-Parametric PV-H (NPPV-H): Non-Parametric Parameter Varying Hammerstein method to identify PV Hammerstein systems.
    %    See E. Sobhani Tehrani and R.E. Kearney, "A Non-Parametric Approach for Identification of Parameter Varying Hammerstein Systems"
    %        IEEE Access, vol. 10, pp. 6348-6362, 11 Jan. 2022, DOI: 10.1109/ACCESS.2022.3141704
    
    properties

    end
    
    methods
        function sys = pvnlbl(a,varargin)
            sys.notes = 'Non-Parametric Parameter Varying Hammerstein (PV-H) model';
            sys.comment = 'Hammerstein cascade elements are modulate by a scheduling variable (SV)';
            
            sys.parameterSet(1)=param('paramName','idMethod', ...
                'paramDefault','',... %'nppv-h'
                'paramHelp', 'This is an iterative method that identifies a non-parameteric model of PV Hammerstein systems',...
                'paramType','select',...
                'paramLimits', {'','nppv-h','sspv-h'});
            
            sys.elements = {pvnl,pvirf};
            
            if nargin==0
                return
            elseif isa(a,'pvnlbl')
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
            methodList= {'','nppv-h','sspv-h'};
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
                case 'nppv-h'
                    j = length(sys.parameterSet);
                    sys.parameterSet(j+1) = param('paramName','irf_len_r',...
                        'paramDefault',0.8,... %0.6,...
                        'paramHelp','Length of IRF in seconds',...
                        'paramType','number',...
                        'paramLimits',[0 inf]);
                    
                    sys.parameterSet(j+2) = param('paramName','nside_r',...
                        'paramDefault',1,...
                        'paramHelp','Number of sides of the IRF',...
                        'paramType','number',...
                        'paramLimits',{1, 2});
                                      
                    %-- Change name from alfa to laguerre_decay
                    sys.parameterSet(j+3) = param('paramName','alfa',...
                        'paramDefault',0.54,... %0.8,...
                        'paramHelp','Constant parameter of the Laguerre expansion',...
                        'paramType','number',...
                        'paramLimits',[0.01 1]); %-- The theoretical lower limit is any number bigger than 0
                    % TBD: Rename parameters to something more easy to
                    % understand
                    %-- The theoretical feasible range for Laguerre exapansion order is [1, inf].
                    %-- But any expansion order < 5 is pratically useless.
                    %-- Any expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from q to laguerre_order
                    sys.parameterSet(j+4) = param('paramName','q',...
                        'paramDefault',8,... %10,...
                        'paramHelp','Order of the Laguerre expansion of the IRF',...
                        'paramType','number',...
                        'paramLimits',[5 16]); 
                                     
                    %-- The theoretical feasible range for SV exapansion order of Laguerre coefficients is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from m to laguerre_sv_order
                    sys.parameterSet(j+5) = param('paramName','m',...
                        'paramDefault',8,... %3,...
                        'paramHelp','Order of SV expansion for IRF Laguerre coefficients',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 corresponds to time-invariant IRF
                    
                    %-- The theoretical feasible range for Input exapansion order of static nl is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from n to nl_input_order
                    sys.parameterSet(j+6) = param('paramName','n',...
                        'paramDefault',9,... %6,...
                        'paramHelp','Order of input expansion for static nonlinearity',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 converts Hammerstein to linear dynamics (TI/PV IRF)
                    
                    %-- The theoretical feasible range for SV exapansion order of static NL coefficients is [0, inf].
                    %-- But expansion order > 16 is computationally too complex and not needed.
                    %-- Change name from p to nl_sv_order
                    sys.parameterSet(j+7) = param('paramName','p',...
                        'paramDefault',7,... %2,...
                        'paramHelp','Order of SV expansion for static nonlinearity coefficients',...
                        'paramType','number',...
                        'paramLimits',[0 16]); %-- Order 0 corresponds to time-invariant IRF
                    
                    sys.parameterSet(j+8) = param('paramName','max_iter',...
                        'paramDefault',500,...
                        'paramHelp','Maximum number of iterations',...
                        'paramType','number',...
                        'paramLimits',[1 1000]);
                    
                    sys.parameterSet(j+9) = param('paramName','threshold',...
                        'paramDefault',10^-10,...
                        'paramHelp','Threshold on SSE for terminitiaon of the iterative search',...
                        'paramType','number',...
                        'paramLimits',[0 1]);
                    
                case 'npnpv-pc'
                    sys.parameterSet = param;
                    j = length(sys.parameterSet);
                    sys.parameterSet(j+1)=param('paramName','idMethod', ...
                        'paramDefault','npnpv-pc',...
                        'paramHelp', 'This is an iterative method that identifies a non-parameteric model of PV Hammerstein systems',...
                        'paramType','select',...
                        'paramLimits', {'npnpv-h','npnpv-pc'});
                    disp('The code for NPNPV-PC method to developed.');
                    
                otherwise
                    %-- If other ident methods are added, put their specific parameters as a new switch case here.    
            end
        end
        
        %% Implementation of various PV Hammerstein identification methods 
        % For example, based on a comment from meeting with Rob, we need something like: nlident(aa,z,'method','npv_hamm')
        % TBD: Make io as z, and sv as a separate input.
        function sys = nlident(sys, z, sv, varargin)
            iMethod = find(strcmp(varargin,'idMethod'));
            newMethod = varargin(iMethod+1);
            %++ ToDo: The above line must be replaced reading the limits from sys object
            %++ something like this:
            % P = sys.parameterSet;
            % methodList = P(1).paramLimits; %-- In principle, we have to read the in index of idMethod instead of assuming it is 1
            methodList= {'nppv-h','sspv-h'};
            if isempty(newMethod) || ~any(strcmp(newMethod, methodList))
                disp('Not a valid identification method option');
                disp (['Available identification method options are: [' strjoin(methodList) ']' ]);
                return
            end
            
            %== Decimation should naturaly be applied to input triplet data (input, sv, output) 
            iDecimation = find(strcmp(varargin,'decimation'));
            decimation = varargin(iDecimation+1); decimation = decimation{1,1};
                       
            switch newMethod{1,1}             
                case 'nppv-h'
                    io = z;      %-- z contains the io data
                    rho = sv;    %-- sv contains the scheduling variable
                    %++ ToDo: The above lines must be removed by updating
                    %         the identification function to accept 
                    %         triplet [input,sv,output] instead of rho as a
                    %         separate signal.
                    
                    %++ ToDo: The following lines shall not be needed 
                    %         when "hammerstein_lpv_laguerre_v02" is
                    %         updated to accept the sys object.
                    %++ For now, extract identification parameters from the sys object outside the function
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
                    static_nl_init = []; 
                    
                    [static_nl,h_r] = hammerstein_lpv_laguerre_v02(io,rho,wNL,wIRF,'irf_len_r',irf_len_r,'n',n,'p',p-1,'m',m,'q',q,'alfa',alfa,'max_iter',max_iter,...
                                                                 'static_nl_param_init',static_nl_param_init,'use_vel_zeroth_exp',use_vel_zeroth_exp,'decimation_ratio',decimation,'static_nl_init',static_nl_init);
                                       
                    %% Setup the objects based on the identification results
                    %++ Extracting SV Tchebychev polynomials from static_nl and setup PVNL 
                    sv_polynoms_pvnl = cell(static_nl.inputExpOrder+1,1);

                    for j = 1:length(sv_polynoms_pvnl)
                        sv_polynoms_pvnl{j,1} = polynom('polyCoef',static_nl.coeffs(j,:)',...
                                                   'polyOrder',static_nl.svExpOrder,...
                                                   'polyType','tcheb');                             
                    end
                    rho_d = decimate_kian(rho,decimation);
                    i_d = decimate_kian(z(:,1),decimation);
                     
                    min_sv = min(rho_d.dataSet); max_sv = max(rho_d.dataSet);
                    min_u = min(i_d.dataSet); max_u = max(i_d.dataSet);
                            
                    PVNL_BASIS = mimobasis;
                    PVNL_BASIS = set(PVNL_BASIS,'coeffs',sv_polynoms_pvnl,...
                                'svRange',[min_sv,max_sv],...
                                'inputRange',[min_u,max_u],...
                                'svExpType','tcheb',...
                                'inputExpType','tcheb',...
                                'svExpOrder',static_nl.svExpOrder,...
                                'inputExpOrder',static_nl.inputExpOrder,...
                                'coeffsStruct',static_nl);
                    
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
                                      'alfa',h_r.LaguerreAlfa,...
                                      'coeffsStruct',h_r);
                    
                    % PVIRF = pvirf; %== TBDel (To Be Deleted) if the following line works
                    PVIRF = pvirf('irfExpansion','laguerre');
                    
                    PVIRF.elements = PVIRF_BASIS;
                    PVIRF.domainIncr = h_r.Ts;
                    
                    %++ Set the cascaded Nonlinearity and Dynamics elements of the identified PV-H system
                    set(sys,'elements',{PVNL,PVIRF});                    
                    
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
                     {'sv_values' nan 'SV values for which to observe the identified system'} ...
            };
            if arg_parse(options,varargin)
                return
            end
            
            %++ Extract elements of the cascade: PVNL and PVIRF
            PVNL = sys.elements{1,1};
            PVIRF = sys.elements{1,2};
            
            if isnan(sv_values)
                subplot(1,2,1)
                plot(PVNL,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv)
                subplot(1,2,2)
                plot(PVIRF,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv)
            else
                subplot(1,2,1)
                plot(PVNL,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv,'sv_values',sv_values)
                subplot(1,2,2)
                plot(PVIRF,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv,'sv_values',sv_values)
            end
        end
        
    end %--> End of methods
end     %--> End of classdef

%% Required functions/subroutines for identification
%- 1. NPNPV-H Identification Method
%- 2. Half wave rectifier fit for initialization of static nonlinearity
%- 3. Data matrix (Gamma_alpha) constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
%- 4. Data matrix (Gamma_c)     constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method

%% 1. NPNPV-H Identification Method
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

%% 2. Half wave rectifier fit for initialization of static nonlinearity
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

%% 3. Data matrix constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
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

%% 4. Data matrix constructor of the Normalized Iterative Convex Search (NICS) algorithm of NPNPV-H identification method
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