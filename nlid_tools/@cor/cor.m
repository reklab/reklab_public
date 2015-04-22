classdef cor < kern
    % cor - correlation class for NLID toolbox
    properties
        Name='cor';
    end
    
    methods
        
        function C = cor (a,varargin)
            %correlation function object
            % Parent: kern
            
            % Copyright 2003, Robert E Kearney and David T Westwick
            % This file is part of the nlid toolbox, and is released under the GNU
            % General Public License For details, see ../copying.txt and ../gpl.txt
            
            j=length(C.parameterSet);
            C.parameterSet(j+1)=param('paramName','biasMode', ...
                'paramDefault', 'unbiased',...
                'paramHelp','Biased or unbiased', ...
                'paramType','select',...
                'paramLimits', {'biased' 'unbiased'});
            C.parameterSet(j+2)=param('paramName','corType',...
                'paramDefault','coeff',...
                'paramHelp','type of correlation function', ...
                'paramType','select',...
                'paramLimits',{'covar' 'correl' 'coeff'});
            
            if nargin==0;
                return
            elseif nargin==1,
                C=nlmkobj(C,a);
            elseif isa(a,'cor')
                C=nlmkobj(a,varargin{:});
            else
               
                C=nlmkobj(C,a,varargin{:});
            end
            
            
        end
        
    end
end

