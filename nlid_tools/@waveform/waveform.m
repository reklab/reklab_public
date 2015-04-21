classdef waveform  < nltop
    % waveform - waveform class for NLID toolbox
    % Generates a variety of waveform  determined by the parameter 'waveformType';
    
    properties
        waveformType=''
        parameterSet=param;
    end
    
    methods
        function sys = waveform  (a,varargin)
            sys.comment='waveform default ';
            if nargin==0;
                sys=nlmkobj(sys,'waveformType','sin');
            elseif nargin==1,
                sys=nlmkobj(sys,a);
            elseif isa(a,'waveform')
                sys=nlmkobj(a,varargin{:});
            else
                sys=nlmkobj(sys,a,varargin{:});
            end
            
        end
        
        function sys= set(sys, varargin)
            sysname=inputname(1);
            % Set method first so that parameters are properly defined for
            % different randvarTypes
            iRVT=find(strcmp(varargin,'waveformType'));
            if ~isempty(iRVT),
                if iRVT<length(varargin),
                    newVarType=varargin{iRVT+1};
                else
                    newVarType='';
                end
                if ~strcmp(sys.waveformType,newVarType),
                    [waveformType,ps]=setwaveformType(newVarType);
                    if isempty(waveformType),
                        return
                    end
                    set@nltop(sys,'waveformType',newVarType);
                    set@nltop(sys,'parameterSet',ps);
                end
            end
            set@nltop(sys, varargin{:});
            assignin('caller',sysname,sys);
        end
        
        
        
        
        
        function Y = nlsim(RV,X)
            % Evalaute a deterministic function at the values of X
            if nargin <2,
                error ('waveform.nlsim requires two inputs: waveform and domain values ');
            end
            if isnumeric(X),
                t=X(:);
            elseif strcmp(class(X),'nldat'),
                t=double(X);
            else
                error ('Second input must be either double or nldat');
            end
            % Check to see if domain is equally sampled
            xRange=range(diff(t));
            if xRange > 10^-6
                domainValues=t;
                domainIncr=NaN;
            else
                domainIncr=mean(diff(t));
                domainValues=NaN;
            end
            Y=nldat;
            set(Y,'comment',['Realization of ' RV.waveformType ' variate']);
            assign(RV.parameterSet);
            switch lower(RV.waveformType)
                case 'impulse'
                    i=min(find (t>= delay));
                    x=zeros(length(t),1);
                    x(i)=1;
                    set(Y,'dataSet',x);
                case 'pulse'
                    iStart=min(find (t>= delay));
                    iEnd = min(find (t>= delay+width));
                    if isempty(iEnd)
                        iEnd=length(t);
                    end
                    x=zeros(length(t),1);
                    for i=iStart:iEnd,
                        x(i)=1;
                    end
                    set(Y,'dataSet',x);
                    
                    
                case 'irf2'
                    a=[ gain damping freq];
                    Y=fgzw(X,a);
                    set(Y,'nLags',length(t), 'irfIdMethod','param');
                    
                case 'ramp'
                    T=1/freq;
                    halfT=T/2.;
                    tMod=mod(t,T);
                    iRise=find(tMod <=halfT);
                    iFall=find(tMod > halfT);
                    x=zeros(length(t),1);
                    x(iRise)=tMod(iRise)*2/T;
                    x(iFall)=1-(tMod(iFall)-halfT)*2/T;
                    set(Y,'dataSet',x);
                case 'sin'
                    x =sin(2*pi*freq*t + pi*phaseShift/180);
                    set(Y,'dataSet',x);
                    
                otherwise
                    error([' waveform type ' detervarType ' not defined']);
            end
            
            set(Y,'domainValues',domainValues,'domainIncr',domainIncr);
        end
        
        function p = nlident(p, z,  varargin)
            % CONSTRUCT an randvar  function object
            
            if nargin==0,
                return
            end
            assign (p.parameterSet);
            tempComment=p.comment;
            if isa(z,'nldat') | isa (z,'double'),
                if isa(z,'double'),
                    z=nldat(z);
                end
                
            end
            
            
            
        end
    end
end
function [waveformType, parameterSet] = setwaveformType (newVarType)
parameterSet=param;
varTypeList={'ramp' 'sin' 'impulse' 'irf2' 'pulse'};
if isempty(newVarType),
    disp('Supported  waveformTypes: ');
    for i=1:length(varTypeList),
        disp(varTypeList{i});
    end
    waveformType='';
    return
elseif ~ismember(lower(newVarType),varTypeList),
    disp([ newVarType ' is not a supported by waveform']);
    disp('Supported  waveformTypes: ');
    for i=1:length(varTypeList),
        disp(varTypeList{i});
    end
    error(['bad waveformType: ' newVarType ]);
    return
end
switch  lower(newVarType)
    case ('impulse')
        waveformType='impulse';
        parameterSet(1)=param('paramName','delay',...
            'paramDefault',0, ...
            'paramHelp','delay)' ,  ...
            'paramType','number',...
            'paramLimits', {-inf inf});
    case ('pulse')
        waveformType='pulse';
        parameterSet(1)=param('paramName','delay',...
            'paramDefault',0, ...
            'paramHelp','delay)' ,  ...
            'paramType','number',...
            'paramLimits', {-inf inf});
        parameterSet(2)=param('paramName','width',...
            'paramDefault',1, ...
            'paramHelp','pulse width )' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
    case{ 'irf2' }
        waveformType='irf2';
        parameterSet(1)=param('paramName','gain',...
            'paramDefault',1, ...
            'paramHelp','gain parameter)' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
        parameterSet(2)=param('paramName','damping',...
            'paramDefault',.5, ...
            'paramHelp','damping parameter)' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
        parameterSet(3)=param('paramName','freq',...
            'paramDefault',2*pi, ...
            'paramHelp','frequency (rad/s)' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
    case ('ramp')
        waveformType='ramp';
        parameterSet(1)=param('paramName','freq',...
            'paramDefault',1, ...
            'paramHelp','frequency (Hz)' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
    case{ 'sin' }
        waveformType='sin';
        parameterSet(1)=param('paramName','freq',...
            'paramDefault',1, ...
            'paramHelp','frequency (Hz)' ,  ...
            'paramType','number',...
            'paramLimits', {0 inf});
        parameterSet(2)=param('paramName','phaseShift',...
            'paramDefault',0, ...
            'paramHelp',' phase (degrees)' ,  ...
            'paramType','number',...
            'paramLimits', {-180 180});
        
        
        
    otherwise
        error(['randvar: unsupported variate type: ' value]);
end
waveformType=newVarType;
comment=[ waveformType ' variate'];
end

