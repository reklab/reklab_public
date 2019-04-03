classdef lnbl < nlm
    % lnbl - linear-nonlinear block model class for NLID toolbox
    
    
    % idMethod - identification method
    %   busgang
    %   hk implements Hunter-Korenberg interation for Wiener cascade models.
    %       See: I.W. Hunter and M.J. Korenberg, The identification of nonlinear 
    %       biological systems: Wiener and Hammerstein cascade models.
    %       Biological Cybernetics, 55:135-144, 1986.
    %   phk implements Paulin-Hunter-Korenberg interation for Wiener cascade models.
% 
%            See M.G. Paulin, A method for constructing data-based models of spiking 
%           neurons using a dynamic linear-static nonlinear cascade.
%           Biological Cybernetics, 69:67-76, 1993.
%   and
%           M.J. Korenberg and I.W. Hunter,  Two methods for identifying Wiener 
%           cascades having noninvertible static nonlinearities.
%           Annals of Biomedical Engineering, 27(6):793-804, 1999.
%   lm -  nonlinear optimization method for Wiener cascade fitting.
% uses a Levenberg-Marquardt second-order gradient descent search.
%
    %
    % initMethod - mehtod used to get inital estimate for linear element 
    %   corel - input-ouput cross-coreelation function
    %   fil - input output impulse response function 
    %   slice2 - random slice of second order correlation function 
    %   slice3 - random slice of third-order correlation function 
    %   eigen - eigen value decomposition of 
    %   gen_eigen - generalize eigne alue decoposition of secondorder
    %   wiener kernel
    
    %% Note: To set paramters of the LE and NLE change them in the elements
    
    properties
        
    end
    
    methods
        function ln = lnbl  (a,varargin)
            ln.parameterSet(1)=param('paramName','idMethod', ...
                'paramDefault','lm',...
                'paramHelp', 'Iterative method used to refine initial estimate',...
                'paramType','select',...
                'paramLimits', {'busgang', 'hk', 'phk', 'lm'});
            setidMethod( ln,'lm');
            ln.comment='LN Model';
            i=irf;
            t=polynom('polyType','tcheb');
            lnps=ln.parameterSet;
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
             
                hlen=get(i,'nLags');
%                 set(i,'nLags',hlen, ...
%                     'irfPseudoInvMode',get(bl,'irfPseudoInvMode'));
                switch P.initMethod
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
                        error(['unrecognized initialization method': initMethod]);
                end
                
                
                x_est = nlsim (h_est,x);
                z1 = cat(2,x_est,y);
                m_est = nlident(p,z1);
                yp = nlsim(m_est,x_est);
                vf = vaf(y,yp);
                set (h_est,'comment','Linear element');
                set (m_est,'comment','Static NonLinear Element');
                set (bl,'elements', { h_est m_est});
                
                
                method=get(bl,'idMethod');
                switch method
                    case 'busgang'
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
        
        function sys= set(sys, varargin)
            sysname=inputname(1);
            % Set method first so that parameters are properly defined
            iMethod=find(strcmp(varargin,'idMethod'));
            if ~isempty(iMethod),
                if length(varargin)>=iMethod+1,
                    sys=setidMethod(sys,varargin{iMethod+1});
                else
                    sys.parameterSet=setval(sys.parameterSet,varargin)
                end
            end
            set@nltop(sys, varargin{:});
            assignin('caller',sysname,sys);
        end
        function sys= setidMethod( ln,newMethod);
            methodList= { 'busgang' 'hk' 'phk' 'lm' };
            if isempty(newMethod) | ~any( strcmp(newMethod, methodList));,
                disp('Not a valid option');
                disp (['Options are: [' strjoin(methodList) ']' ]);
                return
            end
            ln.parameterSet(2)= param('paramName','initMethod', ...
                'paramDefault','eigen',...
                'paramHelp', 'Initial identification method',...
                'paramType','select',...
                'paramLimits',{'correl','fil','slice2','slice3','eigen','gen_eigen'});
            ln.parameterSet(3) = param('paramName','nMaxIts',...
                'paramDefault',20,...
                'paramHelp', 'Maximum number of iterations',...
                'paramType','number',...
                'paramLimits',[1 1000]);
            switch newMethod
                case 'hk'
                    ln.parameterSet(4) = param('paramName','iterationTolerance',...
                        'paramDefault',0.01,...
                        'paramHelp', '% VAF Improvement required to continue HK iteration',...
                        'paramType','number',...
                        'paramLimits', [0 100]);
                    ln.parameterSet=ln.parameterSet(1:4);
                case 'phk'
                    ln.parameterSet(4) = param('paramName','iterationTolerance',...
                        'paramDefault',0.01,...
                        'paramHelp', '% VAF Improvement required to continue HK iteration',...
                        'paramType','number',...
                        'paramLimits', [0 100]);
                    ln.parameterSet(5) = param('paramName','updateGain','paramDefault',1,'paramHelp',...
                        'Gain applied to error in update (alpha)', 'paramType','number',...
                        'paramLimits',[0  inf]);
                    ln.parameterSet=ln.parameterSet(1:5);
                case 'lm'
                    j=length(ln.parameterSet); 
                    ln.parameterSet(j+1)=param('paramName','threshNMSE','paramDefault',.01,'paramHelp',...
                        'NMSE for success','paramType','number','paramLimits',[0 1]);
                    ln.parameterSet(j+2)=param('paramName','accel','paramDefault',.8,'paramHelp',...
                        'ridge multiplied by accell after successful update',...
                        'paramType','number','paramLimits', [0.001 0.999]);
                    ln.parameterSet(j+3)=param('paramName','decel','paramDefault',2,'paramHelp',...
                        'ridge multipled by devel after unsuccessful update',...
                        'paramType','number','paramLimits', [1.0001 inf]);
                    ln.parameterSet(j+4)=param('paramName','delta','paramDefault',10,'paramHelp',...
                        'initial size of ridge added to Hessian','paramType','number',...
                        'paramLimits',[0 inf]);
                otherwise
                    
            end
            
            sys=ln;;
            
            
        end
        
    end
end



