classdef nlmddt < nlm
    % nlmddt - a class for time derivative model (d/dt)
    
    properties
        ddtMethod = '';
    end
    
    methods
        %% nlmddt constructor: Constructs an instance of this class
        function sys = nlmddt(a,varargin)
            sys.comment = 'Numerical time derivative model';
            
            if nargin==0
                return
            elseif nargin==1
                sys = nlmkobj(sys,a);
            elseif isa(a,'nlmddt')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end

            switch sys.ddtMethod
                case ''
                    sys.notes = 'Default numeric method: Fit a parabola with 5-points and use its analytical derivative.';
                otherwise
                    sys.notes = 'The specified numeric method is not implemented yet.'; 
                    disp(sys.notes);
                    return
            end
        end
        
        %% Simulation method for the PVIRF object
        function y = nlsim(sys, u, varargin)
            switch sys.ddtMethod
                case ''
                    y = ddt(u);
                otherwise
                    disp('The specified numeric method is not implemented yet.'); 
                    return
            end
        end
        
    end %--> End of methods
end     %--> End of classdef
