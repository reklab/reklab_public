classdef nlmdelay < nlm
    % nlmddt - a class for time delay model
    
    properties
        delMethod = 'zeropad';
        %-- The nlm class already has a property named inputDelay (in seconds)
    end
    
    methods
        %% nlmdelay constructor: Constructs an instance of this class
        function sys = nlmdelay(a,varargin)
            sys.comment = 'Discrete time delay model';
            
            if nargin==0
                return
            elseif nargin==1
                sys = nlmkobj(sys,a);
            elseif isa(a,'nlmdelay')
                sys = nlmkobj(a,varargin{:});
            else
                sys = nlmkobj(sys,a,varargin{:});
            end

        end
        
        %% Simulation method for the PVIRF object
        function y = nlsim(sys, u, varargin)
            switch sys.delMethod
                case 'zeropad'
                    y = del(u,sys.inputDelay);
                case 'circular'
                    disp(['The ', sys.delMethod, ' method is not implemented yet!'])
                    y = [];
                    return;
                otherwise
                    disp('This delay method does not exist! Available methods are: zeropad and circular')
                    y = [];
                    return
            end
        end
        
    end %--> End of methods
end     %--> End of classdef
