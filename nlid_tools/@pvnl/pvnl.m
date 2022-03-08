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
            };
            if arg_parse(options,varargin)
                return
            end
            
            plot(sys.elements,'n_bins_input',n_bins_input,'n_bins_sv',n_bins_sv);
            xlabel('Input');  
            ylabel('SV'); 
            zlabel('Nonlinearity Output');
            title('PV Static NL');
        end
        
        %% Function to simulate a static PVNL model
        function y = nlsim(sys, z, varargin)
            
        end
        
        %% Function to identify a static PVNL model from data (To be developed)
        function sys = nlident(sys, z, varargin)
             disp('No identification method has yet been developed for PVNL objects');
             return;
        end
        
    end %--> End of methods
end     %--> End of classdef

