classdef F < nldat
    % template for NLID classes
    %   Detailed explanation goes here
    
    properties
        Parameters=param;
    end
    
    methods
        function F = fresp  (a,varargin)
            % Add object  specific parameters
            
            I.Comment='IRF Model';
            if nargin==0;
                return
            elseif nargin==1,
                F=nlmkobj(F,a);
            elseif isa(a,'fresp')
                F=nlmkobj(a,varargin{:});
            else
                F=nlmkobj(f,a,varargin{:});
            end
            
        end
        
      
        
    end
    
end

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt