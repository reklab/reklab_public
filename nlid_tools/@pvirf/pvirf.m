classdef pvirf < pvm
    % pvirf - parameter-varying model class for NLID toolbox, which is 
    % a subclass of the superclass pvm.
    % Note pvirf is a linear parameter varying (LPV) model structure
    
    properties
        domainIncr = 1;
        domainName = 'Lag (s)';
        nSides = 1;
    end
    
    methods
        %% PVIRF Constructor: Constructs an instance of this class
        function sys = pvirf(a,varargin)
            sys.comment = 'PV IRF model';
            
            sys.notes = 'Parameter varying impulse response function (PV IRF)';
            sys.comment = 'Laguerre representation of IRF modulated by a scheduling variable (SV)';
            
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
            mimobasis_irf = set(mimobasis_irf,'svExpType','tcheb','inputExpType','laguerre');
            
            sys.elements = mimobasis_irf;

            if nargin==0
                return
            elseif nargin==1
                sys = nlmkobj(sys,a);
            elseif isa(a,'pvirf')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end
        end
        
        
        %% Function to plot the PV impulse response function (PV IRF)
        function plot(sys,varargin)
            options={{'n_bins_input' 40 'number of bins for a grid on input'} ...
                     {'n_bins_sv' 40 'number of bins for a grid on SV'} ...
            };
            if arg_parse(options,varargin)
                return
            end
            
            plot(sys.elements,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv);
            xlabel('Lags (s)');  
            ylabel('SV'); 
            zlabel('IRF Amplitude');
            title('PV IRF Dynamics');
        end
        
    end %--> End of methods
end     %--> End of classdef

