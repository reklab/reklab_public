classdef  wkern < kern
    % wkern - Wiener kernel class for nLID toolbox
    % Used by wseries.
    
    properties
        name='wkern';
    end
    
    methods
        function WK = wkern (z, varargin)
            iNext=length(WK.parameterSet)+1;
            WK.parameterSet(iNext)=param('paramName','irfMethod',...
                'paramDefault','Toeplitz', ...
                'paramHelp','Identification method' ,...
                'paramType','select',...
                'paramLimits', {'Toeplitz' 'LS'});
            if nargin==0;
                return
            elseif nargin==1,
                WK=nlmkobj(WK,z);
            else

                WK=nlmkobj(WK,z,varargin{:});
            end
            
        end
        
        function y = nlsim ( wk, u, uvar )
            % Simulate response of wkern to input data set
            % input options not fully defined as yet
             order=get(wk,'kernOrder');
            Ts=get(wk,'domainIncr');
            k=get (wk,'dataSet');
            ud=double(u(:,1));
            if nargin < 3
                warning('orthogonalization variance not specified, using signal');
                uvar = std(ud)^2;
            end
            
            % use the Volterra kernel computing routines.
            kv = vkern;
            set(kv,'domainIncr',Ts,'kernOrder',order,'dataSet',k);
            
            
            
            switch order
                case 0
                    y=u*0 + k;
                case 1
                    y = nlsim(kv,u);
                case 2
                    y20 = Ts^2* uvar * sum(diag(k));
                    y = nlsim(kv,u)-y20;
                case 3
                    k31 = sum(squeeze(sum(k)));
                    vk31 = kv;
                    set(vk31,'order',1,'data',k31);
                    y = nlsim(kv,u)-3*uvar*nlsim(vk31,u);
                otherwise
                    error ('nlsim not defined for kernels of order > 3');
            end
            
            
            return
        end
        function i=nlmtst(i)
            %
            % test of wkern identification
            disp('nlmtst for wkern');
            %
            z=rand(10,10);
            i=wkern(z);
            plot (i);
            % wkern/nlmtst
        end
                
        end
    end
    
