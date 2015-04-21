classdef wseries < nlm
    % wseries - Wiener series model class for NLID toolbox
    %  Parent: nlm
    
    % Copyright 1999-2003, Robert E Kearney and David T Westwick
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    
    
    
    properties
    end
    
    methods
        function WS = wseries (z, varargin)
            WS.parameterSet(1) = param('paramName','method', ...
                'paramDefault','LS',...
                'paramHelp','Identifcation Method', ...
                'paramType','select',...
                'paramLimits',{'ls' 'LS' 'Toeplitz' 'toeplitz'});
            
            
            WS.parameterSet(2)=param('paramName','orderMax', ...
                'paramDefault',2,...
                'paramHelp','Maximum order for series' ,...
                'paramType','number', ...
                'paramLimits', {0 3});
            
            WS.parameterSet(3)=param('paramName','nLags', ...
                'paramDefault',NaN,...
                'paramHelp','Number of lags in each kernel' ,...
                'paramType','number',...
                'paramLimits', {0 1000});
            
            WS.parameterSet(4)=param('paramName','variance', ...
                'paramDefault',NaN,...
                'paramHelp','Input variance used in orthogonalization' ,...
                'paramType','number', ...
                'paramLimits', {0 inf});
            
            w0=wkern('kernOrder',0);
            w1=wkern('kernOrder',1);
            w2=wkern('kernOrder',2);
            el = { w0 ; w1; w2 };
            WS.elements=el;
                   
            if nargin==0;
                return
            elseif nargin==1,
                WS=nlmkobj(WS,z);
            else
                args=varargin;
                WS=nlmkobj(WS,z,varargin{:});
            end
        end
        
        
        
        
        
        
    end
    
end

