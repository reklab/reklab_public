classdef nlbl < nlm
    % nlbl - nonlinear-linear block model class for NLID toolbox.
    % Note that the length of the irf and the maximum order of the
    % nonlinearity are parameters of the nl class.
    %
    % idMethods available are: 'hk', 'sls' and 'subspace'
    % HK method is described in:
    %   I.W. Hunter and M.J. Korenberg, The identification of nonlinear
    %   biological systems: Wiener and Hammerstein cascade models.
    %   Biological Cybernetics, 55:135-144, 1986.
    %
    % SLS Method is described in:
    %  See Westwick, D.T., and Kearney, R.E., Separable Least Squares
    %     Identification of Nonlinear Hammerstein Models: Application to
    %     Stretch Reflex Dynamics. Annals of Biomedical Engineering,
    %     29:707--718, 2001.
    %
    % Subspace method is described in:
    %   Jalaleddini, K., and Kearney, R.E., Subspace Identification of SISO Hammerstein Systems:
    %   Application to Stretch Reflex Identification, 2013, in IEEE Transactions on Biomedical Engineering.
    %   Jalaleddini, K., and Kearney R.E., An Iterative Algorithm for the Subspace Identification
    %   SISO Hammerstein Systems, in Proceedings of 18th World Congress IFAC,
    % pp.11779-11784, 2011.
    %
    % Note that the assignment of gain between the nonlinear and noninear
    % is not unique. Methods available to provide consistent results include:
    %   - normCoefNLE which normalizes the magnitdue of the nonlinear elements to one.
    %   - normGainLe which normalizes the steady state gain of the linear element to 1
    %   - normCoefLE. See help on these merthods for details. 
    %
    % Paramters are:
    %   nLagLE - number of lags in the linear element LE, when described by
    %          an IRF
    %   maxOrderNLE - maximum order of polynomial nonlienar element 
    %   threshNSE - threshold of the normalized squared errror for
    %            convergence
    %   displayFlag - display intermediate results
    %
    % Additional parameters for 'sls' method controlling search 
    %    nIterMax - Maximum number of iterations
    %   accel :ridge multiplied by accell after successful update
    %   decel  :ridge multipled by decdel devel after unsuccessful update
    %   delta: initial size of ridge added to Hessian
    %
    % Additional parameters for 'subsapce method' 
     %  hankleSize: Size of hankle matrix (> linear system order)
    %   orderSelectMethodLE :Method for selection of order of linear element
    %       preset - use value of orderLE
    %       manual - interactive manual selection
    %       largest-gap - auotmoatic sletion 
    %   orderLE :Order of linear subsystemn
    %   Delay Input:number of samples of delay between input and output 
    
    
    
       
    % Parent: nlm
    % Parse input
    properties
      idMethod
    end    
    methods
        function nl  = nlbl (z, varargin)
            
            nl.comment='NL Model';
            i=irf;
            t=polynom('polyType','tcheb','polyOrderSelectMode','full');
            nl.elements ={ t i };
           set (nl,'idMethod','hk');
          
            if nargin==0;
                return
            elseif nargin==1,
                nl=nlmkobj(nl,z);
            elseif isa(z,'nlbl')
                nl=nlmkobj(z,varargin{:});
            else
                nl=nlmkobj(nl,z,varargin{:});
            end
            return
        end
        
        function sys= set(sys, varargin)
            sysname=inputname(1); 
            % Set method first so that parameters are proerly defined
            iMethod=find(strcmp(varargin,'idMethod')); 
            if ~isempty(iMethod),
                    [newMethod,ps]=setIdMethod(varargin{iMethod+1});
                    set@nltop(sys,'idMethod',newMethod);
                    set@nltop(sys,'parameterSet',ps);
            end
                    set@nltop(sys, varargin{:});
                      assignin('caller',sysname,sys);
        end
        
        
       
        function nl  = nlident (nl, z, varargin)
            % Identify an NL BLock system
            
            if (nargin > 2),
               set (nl,varargin{:});
            end
            
            if isa(z,'nldat') | isa(z,'double')
                
                if isa(z,'nldat')
                    Ts=z.domainIncr;
                else
                    subsys = nl.elements;
                    f1 = subsys{2};
                    Ts = f1.domainIncr;
                    z = nldat(z,'domainIncr',Ts);
                end
                
                m=get(nl,'idMethod');
                switch lower(m)
                    case 'hk'
                        nl = hammerhk (z,nl,nl.parameterSet);
                        set(nl,'comment',' NL model identified using hk');
                    case 'sls'
                        % Use hk to make initial guess if polynomial has not
                        % been specified
                        % i.e. check to see if its coefficients are NaNs
                        el=nl.elements;
                        coeffs =el{1}.polyCoef;
                        
                        if isnan(coeffs)
                            hk=nlbl;
                            set(hk,'idMethod','hk');
                            % hk=pdefault(hk);                        
                            nl = hammerhk (z,nl, nl.parameterSet);
                        end
                        nl = hammersls(nl,z, nl.parameterSet);
                        set(nl,'comment',' NL model identified using sls');
                    case 'subspace'
                        ztype = class(z);
                        switch ztype
                            case 'segdat'
                                nl = hammer_subspace_short_segment(z,nl.parameterSet);
                            case 'nldat'
                                nl = hammer_subspace(z,nl.parameterSet);
                        end 
                    otherwise
                        disp (['Method ' m ' not defined for nlbl']);
                end
            else
                error('conversions to models of class nlbl not yet implemented');
            end
            nl=normCoefNLE(nl);
        end
        
        function  NN  = normGainLE ( N )
            % Normalize nlbl to assign all gain to the NL element.
            % Scale Hammerstein model apropriately
            % Only valid for low pass systems. 
            
            el=N.elements;
            m=N.elements{1};
            h=N.elements{2};
            hi=cumsum(h)*h.domainIncr;
            gain=double(hi(end));
            if abs(gain)<.01,
                warning('SS Gain is low - may not be highpass');
                gain =0.01*sign(gain);
            end
            m=m.*gain;
            h=h./gain;
            NN=N;
            set(NN,'elements',{ m h});          
        end
        function  NN  = normCoefNLE ( N )
            % Normalize the polynominal coefficeint of a NL model 
            % so that norm(polycoef,2)=1 and 
            %  the first nonsero coefficient  is positive
            m=N.elements{1};
            h=N.elements{2};
            polyCoef=m.polyCoef;
            coefNorm=norm(polyCoef,2);
            iNonZero=find(abs(polyCoef)>.001);
            if ~isempty(iNonZero),
            if polyCoef(iNonZero(1))<0,
                coefNorm=-coefNorm;
            end
            end
            polyCoef=polyCoef/coefNorm;
            set(m,'polyCoef',polyCoef);
            if isa(h,'ssm')
                if ~isnan(coefNorm)
                    set(h,'B',get(h,'B')*coefNorm,'D',get(h,'D')*coefNorm);
                else
                    h = ssm;
                end
            else
                set(h,'dataSet',double(h)*coefNorm);      
            end
            NN=N;
            set(NN,'elements',{ m h});          
        end
        
         function  NN  = normCoefLE ( N )
            % Normalize the IRF coefficeints of a NL model 
            % so that norm(irf,2)=1 and 
            %  the maximum value is positive
            m=N.elements{1};
            h=N.elements{2};
            polyCoef=m.polyCoef;
            if isa(h,'ssm')
                system_ss = ss(h.A,h.B,h.C,h.D,h.domainIncr);
                tf_l = tf(system_ss);
                num = get(tf_l,'num');
                num = num{1};
                den = get(tf_l,'den');
                den = den{1};
                gain = sum(num)/sum(den);
                set(h,'B',get(h,'B')/gain,'D',get(h,'D')/gain);
                set(m,'polyCoef',polyCoef*gain);
            else
                irf=double(h);
                irfNorm=norm(irf,2);
                i=find(abs(irf)==max(abs(irf)));
                if irf(i(1))<0
                    irfNorm=-irfNorm;
                end
                irf=irf/irfNorm;
                polyCoef=polyCoef*irfNorm;
                set(m,'polyCoef',polyCoef);
                set(h,'dataSet',irf); 
            end
            NN=N;
            set(NN,'elements',{ m h});          
        end
        
        function hs = smo (h, npass);
            % smo nlbl object
            npass=1;
            hs=h;
            el=get(h,'elements');
            irf=el{1,2};
            el{1,2}=smo(irf,npass);
            set(hs,'elements',el);
        end
        
        
        
        
    end
end
    
     function [newMethod,ps]= setIdMethod( newMethod);
            methodList= { 'hk' 'sls' 'subspace' };
            if isempty(newMethod) | ~any( strcmp(newMethod, methodList));,
                disp('Not a valid option');
                disp (['Options are: [' strjoin(methodList) ']' ]);
                return
            end                           
            ps=param;
            % Parameters common to all methods         
             ps(1)=param('paramName','nLagLE',...
                        'paramDefault',16,...
                        'paramLimits', [1 inf], ...
                        'paramHelp','Number of lags for linear element IRF ');
                    ps(2)=param('paramName', 'maxOrderNLE',...
                        'paramDefault',12, ...
                         'paramLimits', [1 100], ...
                        'paramHelp','Maximum order of nonlinear element');
                    ps(3)=param('paramName', 'threshNSE',...
                        'paramDefault',.1, ...
                         'paramLimits', [0 inf], ...
                        'paramHelp','threshold on normalized square errors to terminate iteration');
                     ps(4)=param('paramName','displayFlag',...
                        'paramDefault',true,...
                        'paramType','logical', ...
                        'paramHelp','Display flag ');
            % Parameters specific to a method
            switch newMethod
                case 'hk'
                   
                case 'sls'
                    ps(5)=param('paramName','nIterMax',...
                        'paramDefault',20, ...
                         'paramLimits', [1 inf], ...
                        'paramHelp','Maximum number of iterations');
                    ps(6)=param('paramName','accel',...
                        'paramDefault',.8, ...
                         'paramLimits', [0 inf], ...
                        'paramHelp','ridge multiplied by accell after successful update');
                    ps(7)=param('paramName','decel','paramDefault',2, ...
                         'paramLimits', [0 inf], ...
                        'paramHelp','ridge multipled by decdel devel after unsuccessful update');
                    ps(8)=param('paramName','delta',...
                        'paramDefault',10, ...
                         'paramLimits', [0 inf], ...
                        'paramHelp','initial size of ridge added to Hessian');
                case 'subspace'
                     ps(5)=param('paramName','hankleSize',...
                        'paramDefault',20, ...
                         'paramLimits', [1 inf], ...
                        'paramHelp','Size of hankle matrix (> linear system order)');
                    ps(6)=param('paramName','orderSelectMethodLE',...
                        'paramDefault','manual', ...
                        'paramType','select', ...
                        'paramLimits',{'manual' 'largest-gap' 'preset'}, ...
                        'paramHelp','Method for selection of order of linear element');
                    ps(7)=param('paramName','orderLE',...
                        'paramDefault',2, ...
                        'paramType','number', ...
                        'paramLimits',[ 0 inf ], ...
                        'paramHelp','Order of linear subsystem');
                    ps(8)=param('paramName','nDelayInput', ...
                        'paramDefault',0, ...
                         'paramLimits', [1 inf], ...
                        'paramHelp','number of samples of delay ');
                otherwise
                    disp(['specified value of parameter idMethod is not defined:' newMethod]);
                    disp(['Valid options are:' methodList]);
                    error ('Bad option value');
            end
            
        end
        


