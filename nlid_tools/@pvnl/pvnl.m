classdef pvnl < pvm  % pvnl is a subclass of pvm 
    % pvnl - Parameter varying (PV) static nonlinear model class for NLID toolbox.
    
    % pvnl is a nonlinear static model of input and scheduling variable (SV) 
    % in which every coefficient of input nonlinearity (polynomial) is a 
    % static nonlinear (polynomial) function of SV.
    
    % Parent is class pvm 
    
    % Both static nonlinearities (input and SV) are Chebychev polynomials. 
    % Support for additional types of nonlinearities not yet implemented.
    
    % Copyright 2021, Ehsan Sobhani Tehrani and Robert E Kearney
            % This file is part of the nlid toolbox, and is released under the GNU
            % General Public License For details, see ../copying.txt and ../gpl.txt
    
    properties
         
    end
    
    methods
        function sys = pvnl(a,varargin)
            sys.notes = 'Parameter varying static nonlinearity';
            sys.comment = 'Nonlinear function of input modulated by a scheduling variable (SV)';
            
            sys.parameterSet(1) = param('paramName','idMethod', ...
                'paramDefault','',...  % 'paramDefault','npnpv-h'
                'paramHelp', 'This is an iterative method that identifies a non-parameteric model of PV Hammerstein systems',...
                'paramType','select',...
                'paramLimits', {'','npnpv-h','npnpv-pc'});
            
            mimobasis_nl = mimobasis;
            mimobasis_nl = set(mimobasis_nl,'svExpType','tcheb','inputExpType','tcheb');
            
            sys.elements = mimobasis_nl;
            
            if nargin==0
                return
            elseif isa(a,'pvnl')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end
        end
        
        %% Function to plot the PV Nonlinearuty (PV NL)
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
                xlabel('Input');  
                ylabel('SV'); 
                zlabel('Nonlinearity Output');
                title('PV Static NL');
            else
                plot(sys.elements,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv,'sv_values',sv_values);
                xlabel('Input');
                ylabel('Nonlinearity Output');
                nSVs = length(sv_values);
                sv_legends = cell(1,nSVs);
                for i = 1:nSVs
                    sv_legends{1,i} = ['SV=',num2str(sv_values(i))];
                end
                legend(sv_legends);
            end
        end
        
        %% Function to simulate a static PVNL model
        function yp = nlsim(sys, u, sv, varargin)
            PVNL = sys.elements;
            model.static_nl = PVNL.coeffsStruct;
            yp = pvnlSim(model,u,sv);
        end
        
        %% Function to get the snapshots of the PVNL model
        function snapshotModels = snapshot(sys,sv_values)
            mimoSnapshots = snapshot(sys.elements,sv_values);
            snapshotModels.svValues = sv_values;
            nModels = length(sv_values);
            snapshotModels.elements = cell(1,nModels);
            for j = 1:nModels
                x = nldat(mimoSnapshots.domainValues);
                y = nldat(mimoSnapshots.dataSet(:,j));
                z = cat(2,x,y);
                element = polynom('polyType','tcheb','polyOrder',sys.elements.inputExpOrder);
                element = nlident(element,z);
                snapshotModels.elements{1,j} = element;
            end
            
        end
        
        %% Function to identify a static PVNL model from data (To be developed)
        function sys = nlident(sys, z, varargin)
             disp('No identification method has yet been developed for PVNL objects');
             return;
        end
        
    end %--> End of methods
end     %--> End of classdef

%% 1. PV-NL simulation function
function yp = pvnlSim(model,u,rho)

%=== Reading data parameters
ts = u.domainIncr;
fs = 1/ts;
nSamples = length(u.dataSet);

%=== Reading model parameters
static_nl = model.static_nl;

avg_i = static_nl.inputNormalization(1);
rng_i = static_nl.inputNormalization(2);

avg_rho = static_nl.svNormalization(1);
rng_rho = static_nl.svNormalization(2);

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

yp = nldat(z,'domainIncr',ts);

end