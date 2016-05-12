classdef vseries < nlm
    % vseries - Volterra series model class for NLID toolbox
    %  Parent: nlm
    
    % Copyright 1999-2003, Robert E Kearney and David T Westwick
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    properties
    end
    
    methods
        function VS = vseries (z, varargin)
            % Volterra Series Model
            %  Parent: nlm
            
            VS.parameterSet(1) = param('paramName','idMethod', ...
                'paramDefault','foa',...
                'paramHelp','Identification method' , ...
                'paramType','select',...
                'paramLimits', {'FOA' 'foa' 'Laguerre' 'laguerre'});
            VS.parameterSet(2)=param('paramName','nLags', ...
                'paramDefault',16,...
                'paramHelp','Number of lags in each kernel' ,...
                'paramType','number', ...
                'paramLimits', [0 1000]);
            VS.parameterSet(3)=param('paramName','vsOrderMax', ...
                'paramDefault',2, ...
                'paramHelp','Maximum order for series' ,...
                'paramType','number', ...
                'paramLimits', [0 2]);
            
            VS.parameterSet(4) = param('paramName','alphaLaguerre', ...
                'paramDefault',NaN,...
                'paramHelp','decay parameter for Laguerre filters',...
                'paramType','number', ...
                'paramLimits',[eps 1-eps]);
            VS.parameterSet(5) = param('paramName','nLaguerreFilt',...
                'paramDefault',NaN,...
                'paramHelp','Number of Laguerre filters in expansion',...
                'paramType','number',...
                'paramLimits',[1 100]);
            
            VS.parameterSet(6) = param('paramName','delayLaguerre' ,...
                'paramDefault',0,...
                'paramHelp','delay before start of Laguerre fitlers',...
                'paramType','number',...
                'paramLimits',[0 1000]);
            
            % Initialize elements
            v0=vkern('kernOrder',0);
            v1=vkern('kernOrder',1);
            v2=vkern('kernOrder',2);
            el = { v0 ; v1; v2 };
            
            set(VS,'elements',el);
            
            if nargin==0;
                return
            elseif nargin==1,
                VS=nlmkobj(VS,z);
            else
                args=varargin;
                VS=nlmkobj(VS,z,varargin{:});
            end
        end
        
        function vs = nlident (vs, z, varargin)
            % Identify a Volterra series model
            
            if nargin > 2,
                set(vs,varargin{:});
            end
            
            if isa(z,'nldat') | isa(z,'double')
                
                if isa(z,'nldat')
                    Ts=get(z,'domainIncr');
                else
                    stuff = get(vs,'elements');
                    vk0 = stuff{1};
                    Ts = get(vk0,'domaInincr');
                end
                
                zd=double(z);
                u=zd(:,1);
                y=zd(:,2);
                numlags=get(vs,'nLags');
                if isnan(numlags),
                    numlags=min (16,length(z)/100);
                    set(vs,'nLags',numlags);
                end
                %
                % Compute first three kernels uisng fast orthongal method
                method = get(vs,'idMethod');
                switch lower(method)
                    case 'foa'
                        vs = foa(vs,u,y,Ts);
                    case 'laguerre'
                        vs = lag_id(vs,u,y,Ts);
                end
            elseif isa(z,'nlbl'),
                vs=nl2vs(z,vs);
            elseif isa(z,'lnbl'),
                vs=ln2vs(z,vs);
            elseif isa(z,'lnlbl'),
                vs=lnl2vs(z,vs);
            elseif isa(z,'wseries')
                vs = ws2vs(z,vs);
            elseif isa(z,'wbose')
                vs = wb2vs(z,vs);
            elseif isa(z,'pcascade'),
                vs=pc2vs(z,vs);
            elseif isa(z,'nlm'),
                vs = nlm2vs(z,vs);
            else
                error (['vseries does not support inputs of type:' class(z)] );
            end
            
            function vs = pc2vs (pc, vsin);
                % Generate Volterra series from parallel cascade description
                %
                ordermax=get(vsin,'vsOrderMax');
                k=cell(ordermax+1,1);
                for order=0:ordermax
                    vskern = pcas2volt(pc,order);
                    k{order+1,1}=vskern;
                end
                
                vs=vsin;
                set(vs,'elements',k);
            end
        end
        
        
        function y = nlsim ( vs, u )
            % Simulate response of a Volterra series model  to an input signal
            % input options not fully defined as yet
            
            %Q =get (vs,'order');
            kernels = get(vs,'elements');
            num_kerns = length(kernels);
            ud=double(u);
            y=nldat(u);
            yd = zeros(size(ud));
            
            for i = 1:num_kerns
                vk = kernels{i};
                yd = yd + double(nlsim(vk,u));
            end
            set(y,'dataSet',yd);
        end
        
        
        
        function plot (vs)
            % Plot a volterra series model
            %
            kernels =get(vs,'elements');
            Qmax = length(kernels)-1;
            
            switch Qmax
                case 0
                    subplot(111);
                    plot(kernels{1});
                case 1
                    subplot(211)
                    plot(kernels{1});
                    subplot(212)
                    plot(kernels{2});
                case 2
                    subplot(221)
                    plot(kernels{1});
                    subplot(222)
                    plot(kernels{2});
                    subplot(212)
                    plot(kernels{3});
                otherwise
                    subplot(221)
                    plot(kernels{1});
                    subplot(222)
                    plot(kernels{2});
                    subplot(223)
                    plot(kernels{3});
                    subplot(224)
                    plot(kernels{4});
            end
            
            streamer (get(vs,'comment'));
        end
    end
    
    
    
    
end
