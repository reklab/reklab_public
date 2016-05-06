classdef randvar  < nltop
    % randvar - random variate class for NLID toolbox
    % Generates variates from a variety of random distirubtions defined by
    % paramter 'randVarType';
    
    properties
        randvarType=''
        parameterSet=param;
    end
    
    methods
        function sys = randvar  (a,varargin)
            sys.comment='randvar default ';
            if nargin==0;
                sys=nlmkobj(sys,'randvarType','Normal');
            elseif nargin==1,
                sys=nlmkobj(sys,a);
            elseif isa(a,'randvar')
                sys=nlmkobj(a,varargin{:});
            else
                sys=nlmkobj(sys,a,varargin{:});
            end
            
        end
        
        function sys= set(sys, varargin)
            sysname=inputname(1);
            % Set method first so that parameters are properly defined for
            % different randvarTypes
            iRVT=find(strcmp(varargin,'randvarType'));
            if ~isempty(iRVT),
                if iRVT<length(varargin),
                    newVarType=varargin{iRVT+1};
                else
                    newVarType='';
                end
                if ~strcmp(sys.randvarType,newVarType),
                    [randVarType,ps]=setRandVarType(newVarType);
                    if isempty(randVarType),
                        return
                    end
                    set@nltop(sys,'randvarType',newVarType);
                    set@nltop(sys,'parameterSet',ps);
                end
            end
            set@nltop(sys, varargin{:});
            assignin('caller',sysname,sys);
        end
        
        function x = nlident( r,z,varargin)
            if nargin < 2,
                disp('nlident for ransdvar takes two inputs : radnvar, domain' );
            elseif nargin > 2,
                set(i,varargin{:});
            end
            x=nlsim(r,z);
        end
        
        
        
        function Y =nlsim(RV,X)
            % Generate a realization of a randvar with the same dimensions
            % as X
            if nargin <2,
                error ('randvar.nlsim requires two inputs randvar and nldat ');
            end
            if isnumeric(X),
                X=nldat(X);
            elseif ~strcmp(class(X),'nldat'),
                error ('Second input must be either double or nldat');
            end
            [m,n,p]=size(X);
            Y=X;
            set(Y,'comment',['Realization of ' RV.randvarType ' variate']);
            pList= cat(2, getParamValCell(RV.parameterSet), m,n,p);
            set(Y,'dataSet',random(RV.randvarType,pList{:}));
        end
        
        
        
        
    end
end

function [randvarType, parameterSet] = setRandVarType (newVarType)
parameterSet=param;
varTypeList={'Binomial'  'Chisquare' 'Exponential' 'F' 'Normal'    'Lognormal' 'Poisson' 'T' 'Uniform' 'Discrete Uniform' };
if isempty(newVarType),
    disp('Supported  randvarTypes: ');
    for i=1:length(varTypeList),
        disp(varTypeList{i});
    end
    randvarType='';
    return
elseif ~ismember(newVarType,varTypeList),
    disp([ newVarType ' is not a supported by randvar']);
    disp('Supported  randvarTypes: ');
    for i=1:length(varTypeList),
        disp(varTypeList{i});
    end
    error('bad randvarType');
    return
end
switch  lower(newVarType)
    case ('binomial')
        randvarType='Binomial';
        parameterSet(1)=param('paramName','nTrials',...
            'paramDefault',10, ...
            'paramHelp','' ,  ...
            'paramType','number',...
            'paramLimits', [1 inf]);
        parameterSet(2)=param('paramName','prob',...
            'paramDefault',.5, ...
            'paramHelp','Proabilty of Success' ,  ...
            'paramType','number',...
            'paramLimits', [0 1.]);
    case{ 'chisquare', 'chi2' }
        randvarType='Chisquare';
        parameterSet(1)=param('paramName','dof',...
            'paramDefault',5, ...
            'paramHelp','degrees of freedom' ,  ...
            'paramType','number',...
            'paramLimits', [0 inf]);
    case {'exponential' , 'exp'}
        randvarType='Exponential';
        parameterSet(1)=param('paramName','mean',...
            'paramDefault',1, ...
            'paramHelp','mean value' ,  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
    case {'F' , 'f'}
        randvarType='F';
        parameterSet(1)=param('paramName','DOF1',...
            'paramDefault',5, ...
            'paramHelp','First degree of freedom' ,  ...
            'paramType','number',...
            'paramLimits', [1 20]);
        parameterSet(2)=param('paramName','DOF2',...
            'paramDefault',3, ...
            'paramHelp','Second degree of freedom' ,  ...
            'paramType','number',...
            'paramLimits', [1 20]);
    case{ 'normal'  , 'norm' }
        randvarType='Normal';
        parameterSet(1)=param('paramName','mean',...
            'paramDefault',0, ...
            'paramHelp','mean value' ,  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
        parameterSet(2)=param('paramName','std',...
            'paramDefault',1, ...
            'paramHelp','Standard Deviation',  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
    case{ 'lognormal'  , 'logn' }
        randvarType='Lognormal';
        parameterSet(1)=param('paramName','mean',...
            'paramDefault',0, ...
            'paramHelp','mean value' ,  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
        parameterSet(2)=param('paramName','std',...
            'paramDefault',1, ...
            'paramHelp','Standard Deviation',  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
    case{ 'poisson', 'pois' }
        randvarType='Poisson';
        parameterSet(1)=param('paramName','mean',...
            'paramDefault',1, ...
            'paramHelp','mean value' ,  ...
            'paramType','number',...
            'paramLimits', [0 inf]);
    case{ 't' }
        randvarType='T';
        parameterSet(1)=param('paramName','dof',...
            'paramDefault',10, ...
            'paramHelp','mean value' ,  ...
            'paramType','number',...
            'paramLimits', [0 inf]);
    case {'uniform' , 'uni'}
        randvarType='uniform';
        parameterSet(1)=param('paramName','minVal',...
            'paramDefault',0, ...
            'paramHelp','min value' ,  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
        parameterSet(2)=param('paramName','maxVal',...
            'paramDefault',1, ...
            'paramHelp','maximum value',  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
        
    case {'discrete uniform' , 'unid'}
        randvarType='Discrete Uniform';
        parameterSet(1)=param('paramName','maxVal',...
            'paramDefault',10, ...
            'paramHelp','min value' ,  ...
            'paramType','number',...
            'paramLimits', [-inf inf]);
        
    otherwise
        varTypeList={ 'Chisquare' 'Exponential' 'Normal'    'Lognormal' 'Poisson' 'T' 'Uniform' 'Discrete Uniform' }
        error(['randvar: unsupported variate type: ' value]);
end
comment=[ randvarType ' variate'];
end

