classdef lnbl < nlm
    % lnbl - linear-nonlinear block model class for NLID toolbox
    
    
    properties
        
    end
    
    methods
        function ln = lnbl  (a,varargin)
            % Add object  specific parameters
            ln.parameterSet(1)=param('paramName','lnIdMethod', ...
                'paramDefault','lm',...
                'paramHelp', 'Iterative method used to refine initial estimate',...
                'paramType','select',...
                'paramLimits', {'bussgang', 'hk', 'phk', 'lm'});
            ln.parameterSet(2)= param('paramName','lnInitMethod', ...
                'paramDefault','eigen',...
                'paramHelp', 'Initial identification method',...
                'paramType','select',...
                'paramLimits',{'correl','fil','slice2','slice3','eigen','gen_eigen'});
            ln.parameterSet(3) = param('paramName','nMaxIts',...
                'paramDefault',20,...
                'paramHelp', 'Maximum number of iterations',...
                'paramType','number',...
                'paramLimits',[1 1000]);
            ln.parameterSet(4) = param('paramName','iterationTolerance',...
                'paramDefault',0.01,...
                'paramHelp', '% VAF Improvement required to continue HK iteration',...
                'paramType','number',...
                'paramLimits', [0 100]);
            
            ln.parameterSet(5) = param('paramName','updateGain','paramDefault',1,'paramHelp',...
                'Gain applied to error in update (alpha)', 'paramType','number',...
                'paramLimits',[0  ]});
            
            %  lm parameters
            
            ln.parameterSet(6)=param('paramName','threshNMSE','paramDefault',.01,'paramHelp',...
                'NMSE for success','paramType','number','paramLimits',[0 1]);
            ln.parameterSet(7)=param('paramName','accel','paramDefault',.8,'paramHelp',...
                'ridge multiplied by accell after successful update',...
                'paramType','number','paramLimits', [0.001 0.999]);
            ln.parameterSet(8)=param('paramName','decel','paramDefault',2,'paramHelp',...
                'ridge multipled by devel after unsuccessful update',...
                'paramType','number','paramLimits', [1.0001 inf]);
            ln.parameterSet(9)=param('paramName','delta','paramDefault',10,'paramHelp',...
                'initial size of ridge added to Hessian','paramType','number',...
                'paramLimits',[0 inf]);
            ln.comment='LN Model';
            i=irf;
            t=polynom('polyType','tcheb');
            lnps=ln.parameterSet;
            ips=get(i,'parameterSet');
            pps=get(t,'parameterSet');
            L=[ lnps ips pps];
            ln.parameterSet=L;
            set(ln,'elements', { i t });
            if nargin==0;
                return
            elseif nargin==1,
                ln=nlmkobj(ln,a);
            elseif isa(a,'lnbl')
                ln=nlmkobj(a,varargin{:});
            else
                ln=nlmkobj(ln ,a,varargin{:});
            end
            
        end
        
        function bl  = nlident (bl, z, varargin)
            % Identify a lnbl
            if nargin > 2
                set (bl,varargin{:});
            end
            if isa(z,'nldat') | isa(z,'double')
                if isa(z,'nldat')
                    Ts=z.domainIncr;
                else
                    subsys = get(bl,'elements');
                    f1 = subsys{1};
                    Ts = f1.domainIncr;
                    z = nldat(z,'domainIncr',Ts);
                end
                subsys = bl.elements;
                P=getParamValStruct(bl.parameterSet); 
                i = subsys{1,1};  % IRF
                p = subsys{1,2};  % Polynomial
                x=z(:,1);
                y=z(:,2);
                hlen=get(bl,'nLags');
                set(i,'nLags',hlen, ...
                    'irfPseudoInvMode',get(bl,'irfPseudoInvMode'));
                switch P.lnInitMethod
                    case 'correl'
                        phi = cor('kernOrder',1,'nLags',hlen);
                        phi = nlident(phi,z);
                        h_est = i;
                        set(h_est,'dataSet',double(phi));
                    case 'fil'
                        h_est = nlident(i, z);
                    case 'slice2'
                        h_est = cor_slice(i,z,2);
                    case 'slice3'
                        h_est = cor_slice(i,z,3);
                    case 'eigen'
                        phi = cor;
                        set(phi,'kernOrder',2,'nLags',hlen);
                        phi = nlident(phi,z);
                        [U,S,V] = svd(double(phi));
                        h_est = i;
                        set(h_est,'dataSet',U(:,1),'domainIncr',z.domainIncr);
                    case 'gen_eig'
                        h_est = wiener_2(z,i);
                    otherwise
                        error(['unrecognized initialization method': lnInitMethod]);
                end
                
                
                x_est = nlsim (h_est,x);
                z1 = cat(2,x_est,y);
                set (p,'polyOrderMax',get(bl,'polyOrderMax'),...
                    'polyOrderSelectMode', get(bl,'polyOrderSelectMode'));
                m_est = nlident(p,z1);
                yp = nlsim(m_est,x_est);
                vf = vaf(y,yp);
                set (h_est,'comment','Linear element');
                set (m_est,'comment','Static NonLinear Element');
                set (bl,'elements', { h_est m_est});
                
                
                method=get(bl,'lnIdMethod');
                switch method
                    case 'bussgang'
                        set(bl,'comment','LN model identified using Busgang''s theorm');
                    case 'hk'
                        bl = hk_ident(bl,z);
                    case 'phk'
                        bl = phk_ident(bl,z);
                    case 'lm'
                        bl = lm_ident(bl,z);
                    otherwise
                        error(['unrecognized identification method:' method]);
                end
                
            else
                error('conversions to models of class lnbl not yet implemented');
            end
            function h = cor_slice(h,z,order)
                
              P=getParamValStruct(h.parameterSet);
                ud = double(z(:,1));
                yd = double(z(:,2));
                N = length(ud);
                uny = z;
                hlen=P.nLags;
                
                switch order
                    case 2
                        lag = floor(hlen*rand(1));
                        udel = [zeros(lag,1);ud(1:N-lag)];
                        set(uny,'dataSet',[ud udel.*yd]);
                    case 3
                        lag1 = floor(hlen*rand(1));
                        ud1 = [zeros(lag1,1);ud(1:N-lag1)];
                        lag2 = floor(hlen*rand(1));
                        ud2 = [zeros(lag2,1);ud(1:N-lag2)];
                        set(uny,'dataSet',[ud ud1.*ud1.*yd]);
                    otherwise
                        error('unsupported slice order');
                end
                
                phi = cor('kernOrder',1,'nLags',hlen);
                phi = nlident(phi,uny);
                hd = phi.dataSet;
                
                
                gain = std(yd);
                switch order
                    case 2
                        hd(lag+1) = hd(lag+1) + randn(1)*gain;
                    case 3
                        hd(lag1+1) = hd(lag1+1) + randn(1)*gain;
                        hd(lag2+1) = hd(lag2+1) + randn(1)*gain;
                    otherwise
                        error('unsupported slice order');
                end
                set(h,'dataSet',hd);
            end
        end
    end
    
end


