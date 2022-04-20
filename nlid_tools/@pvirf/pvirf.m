classdef pvirf < pvm
    % pvirf - parameter-varying model class for NLID toolbox, which is 
    % a subclass of the superclass pvm.
    % Note pvirf is a linear parameter varying (LPV) model structure
    
    properties
        domainIncr = 1;
        domainName = 'Lag (s)';
        nSides = 1;
        irfExpansion = '';
    end
    
    methods
        %% PVIRF Constructor: Constructs an instance of this class
        function sys = pvirf(a,varargin)
            sys.comment = 'PV IRF model';
            sys.notes = 'Parameter varying impulse response function (PV IRF)';
            
            if nargin==0
                return
            elseif nargin==1
                sys = nlmkobj(sys,a);
            elseif isa(a,'pvirf')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end

            switch sys.irfExpansion
                case 'laguerre'
                    sys.comment = 'Laguerre representation of IRF modulated by a scheduling variable (SV)';
                case 'none'
                    sys.comment = 'IRF coefficients modulated by a scheduling variable (SV)';
                otherwise
                    sys.comment = 'The specified IRF expansion is not supported. Supported types are: laguerre and none (i.e., use IRF coefficients)'; 
                    disp(sys.comment);
                    return
            end
            
            % add object specific parameters
%             sys.parameterSet(1) =param('paramName','irfBasisExpansion','paramDefault','laguerre', ...
%                 'paramHelp','representation of IRF using coefficients or with basis expansions',...
%                 'paramType','select','paramLimits',{'none','tcheb'});
%             sys.parameterSet(2) =param('paramName','svPolyType','paramDefault','tcheb', ...
%                 'paramHelp','sv polynomial expansion type',...
%                 'paramType','select','paramLimits',{'tcheb'});
%             sys.parameterSet(3) =param('paramName','orderSelectMode','paramDefault','manual',...
%                 'paramHelp','order selection method for expansions of IRF basis and SV expansion',...
%                 'paramType','select','paramLimits',{'manual','auto'});
            
            mimobasis_irf = mimobasis;
            
            switch sys.irfExpansion
                case 'laguerre'
                    mimobasis_irf = set(mimobasis_irf,'svExpType','tcheb','inputExpType','laguerre');
                case 'none'
                    mimobasis_irf = set(mimobasis_irf,'svExpType','tcheb','inputExpType','none');
                    mimobasis_irf = set(mimobasis_irf,'nSides',sys.nSides);
                otherwise
                    disp('This IRF expansion type is not supported. Supported types are: laguerre and none (i.e., uses IRF coefficients).')
                    return
            end    
            
            sys.elements = mimobasis_irf;

        end
        
        
        %% Function to plot the PV impulse response function (PV IRF)
        function plot(sys,varargin)
            options={{'n_bins_input' 40 'number of bins for a grid on input'} ...
                     {'n_bins_sv' 40 'number of bins for a grid on SV'} ...
                     {'sv_values' nan 'SV values for which to observe the identified system'} ...
            };
        
            if arg_parse(options,varargin)
                return
            end
            
            if isnan(sv_values)
                plot(sys.elements,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv);
                xlabel('Lags (s)');  
                ylabel('SV'); 
                zlabel('IRF Amplitude');
                title('PV IRF Dynamics');
            else
                plot(sys.elements,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv,'sv_values',sv_values);
                xlabel('Lags (s)')
                ylabel('IRF Amplitude')
                nSVs = length(sv_values);
                sv_legends = cell(1,nSVs);
                for i = 1:nSVs
                    sv_legends{1,i} = ['SV=',num2str(sv_values(i))];
                end
                legend(sv_legends);
            end
        end
        
        %% Simulation method for the PVIRF object
        function yp = nlsim(sys, u, sv, varargin)
            PVIRF = sys.elements;
            %++ PVIRF output simulation depends on the IRF expnasion type
            switch sys.irfExpansion              
                case 'laguerre'
                    model.h_r = PVIRF.coeffsStruct;
                    yp = pvlaguerreSim(model,u,sv);
                
                case 'none'
                    model.h_i = PVIRF.coeffsStruct;
                    yp = pvirfSim(model,u,sv);
                    
                otherwise
                    disp('This input expansion type is not supported. Supported types are: laguerre and none (i.e., IRF coefficients).')
                    yp = NaN;
                    return
            end
        end
        
        %% Function to get the snapshots of the PVIRF model
        function snapshotModels = snapshot(sys,sv_values)
            mimoSnapshots = snapshot(sys.elements,sv_values);
            snapshotModels.svValues = sv_values;
            nModels = length(sv_values);
            snapshotModels.elements = cell(1,nModels);
            for j = 1:nModels
                element = irf('domainIncr',sys.domainIncr,...
                              'domainValues',mimoSnapshots.domainValues,...
                              'dataSet',mimoSnapshots.dataSet(:,j));
                snapshotModels.elements{1,j} = element;
            end
            
        end
        
    end %--> End of methods
end     %--> End of classdef

%% 1. PV-Laguerre simulation function
function yp = pvlaguerreSim(model,u,rho)

%=== assign u as nldat object to z as array (TBR = To Be Replaced)
nSamples = length(u.dataSet);
z = u.dataSet;

%=== Reading model parameters
h_r = model.h_r;

%=== Reading data parameters
ts = h_r.Ts;
fs = 1/ts;
if ts ~= u.domainIncr
    disp('The sampling time of the input does not match that of the IRF dynamics of the Hammerstein cascade!')
    yp = [];
    return;
end

%=== Normalizing the SV
avg_rho = h_r.svNormalization(1);
rng_rho = h_r.svNormalization(2);

rho_n = (rho.dataSet - avg_rho)*2/rng_rho;

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

%% 2. PV-IRF simulation function
function yp = pvirfSim(model,u,rho)

%=== Reading data parameters
ts = u.domainIncr;
fs = 1/ts;
nSamples = length(u.dataSet);

%=== Reading model parameters
h_i = model.h_i;
intrinsic = h_i.data;
nside_i = h_i.nSides;
nLags_i = h_i.nLags;
p_i = h_i.svExpOrder;

%=== Normalizing the SV
avg_rho = h_i.normalization(1);
rng_rho = h_i.normalization(2);

rho_n = (rho.dataSet - avg_rho)*2/rng_rho;

%=== Computing the lagged positions for the Intrinsic IRF 
Li = (nLags_i - 1) / nside_i;

if nside_i == 1
    lags_i = ts * (0:Li);
else
    lags_i = ts * (-Li:1:Li);
end

u_lags = zeros(nSamples,length(nLags_i));
for d = 1:nLags_i
    u_lags(:,d) = del(u.dataSet,fs,lags_i(d));
end

%== Computing the Tchebychev expansions of scheduling for intrinsic pathway
Gpi =  multi_tcheb(rho_n,p_i);

%== The Regressor of the Intrinsic Pathway
Phi_intrinsic_expand = zeros(nSamples,(p_i+1)*(nLags_i));

col = 0;
for d = 1:nLags_i
    for pic = 1:p_i+1
        col = col + 1;
        Phi_intrinsic_expand(:,col) = u_lags(:,d).*Gpi(:,pic);        
    end
end

%== Intrinsic Torque
yp = Phi_intrinsic_expand*intrinsic;
yp = nldat(yp,'domainIncr',ts);

end