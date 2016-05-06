classdef param
    % param - parameter class for NLID toolbox
    
    % param objects use to hold user-defined properties of an object than
    % may vary from object to object. 

    
    properties
        paramDefault= NaN  ;
        paramHelp='Help';
        paramLimits=[-inf inf]; % Specify as a vecctor  if paramtype is number or string array if select
        paramName='Param Name';
        paramType='number'; 
        paramValue=NaN;
    end
    
    methods
        function P = param(a,varargin)
            if nargin==0;
                return
            elseif nargin==1,
                P=nlmkobj(P,a);
            else
                args=cat(2,{a},varargin);
                set(P,args{:});
            end
        end
        
        function out = set (P, varargin)
            % Set a parameter value
            name = inputname(1);
            nIn=nargin-1;
            if isempty(name),
                error('First argument to SET must be a named variable.')
            elseif rem(nIn,2)~=0,
              disp(P)
            end
            
            for i=1:2:nIn,
                Property = varargin{i};
                if ~any(strcmp(Property, fieldnames(P))),
                    error (['Not a valid parameter field: ' Property]);
                end
                Value = varargin{i+1};
                if strcmp(Property,'paramType'),
                    % check for valid types
                    validTypes= { 'number' 'string' 'logical' 'select' };
                    if ~any(strcmp(Value, validTypes)),
                        disp(['Valid paramTypes:' validTypes]);
                        error ([Value ' is not a valid paramType' ]);
                    end
                end
                P.(Property) = Value;
            end % for
            assignin('caller',name,P)
            
        end
        
        function assign (p)
            % assign values of parameters in calling workspace
            n=length(p);
            for i=1:length(p)
                name=p(i).paramName;
                val=value(p(i));
                assignin('caller',name,val)
            end
        end
        
        function disp(P),
            nParam= length(P);
            for i=1:nParam,
                pVal=getParamValue(P,P(i).paramName);
                     disp(['  ' P(i).paramName ':' num2str(pVal(:)')]); 

            end
        end
       
        
        function dispFull(P)
            l=length(P);
            pad=blanks(4);
            for i=1:l
                disp(['paramName: ' P(i).paramName]);               
                disp(['paramHelp:' P(i).paramHelp]);
                disp(['paramType:' P(i).paramType]);
                limits = P(i).paramLimits;
                ltype = P(i).paramType;
                switch lower(ltype)
                    case 'number'
                         disp(['paramValue:' num2str(P(i).paramValue)]); 
                          disp(['paramDefault:' num2str(P(i).paramDefault)]); 
                        llim = num2str(limits(1));
                        ulim = num2str(limits(2));
                        disp(['paramLimits: [ ' llim pad ulim ' ]']);
                    case 'select'
                        val=(P(i).paramValue);
                        if isnan(val),
                            val='NaN';
                        end
                         disp(['paramValue:' val]); 
                         disp(['paramDefault:'  P(i).paramDefault]), 
                        nopt = length(limits);
                        lstr = ' ';
                        for j = 1:nopt-1
                            lstr = [lstr limits{j} ', '];
                        end
                        lstr = [lstr limits{nopt}];
                        disp(['paramLimits: [ ' lstr ' ]']);
                    case 'string'
                         disp(['paramValue:' (P(i).paramValue{1})]); 
                end 
                 disp('-');
            end
           
        end
        
        function Pout = getParamValue ( P, paramName)
            j=pindex(P,paramName);
            if j==0
                warning (['Parameter name does not exist:' paramName]);
                Pout=NaN;
            else
                if isnan(P(j).paramValue),
                    Pout=P(j).paramDefault;
                else
                    Pout=P(j).paramValue;
                end
            end
        end
        
        function Pout = getParamValStruct ( P )
            
            for i=1:length(P),
                outName=P(i).paramName;
                
                if isnan(P(i).paramValue{1}),
                    Pout.(outName)=P(i).paramDefault;
                else
                    Pout.(outName)=P(i).paramValue{:};
                end
            end
        end
        
        function Pout = getParamValCell ( P )
            % Return a cell array of parameter values
            for i=1:length(P),
                if isnan(P(i).paramValue),
                    Pout{i}=P(i).paramDefault;
                else
                    Pout{i}=P(i).paramValue;
                end
            end
        end
        
          function pOut = getParamArgList ( P )
            % Return a cell array of parameter values for use as an
            % argument list 
            j=0;
            for i=1:length(P),
                j=j+1;
                pOut{j}=P(i).paramName;
                j=j+1;
                if isnan(P(i).paramValue{1}),
                    pOut{j}=P(i).paramDefault;
                else
                    pOut{j}=P(i).paramValue{:};
                end
            end
        end
        
        function j = pindex (P, S);
            % Determine the index for a parameter name
            j=0;
            Pout=param;
            for i=1:length(P),
                if strcmp(P(i).paramName,S)
                    j=i;
                end
            end
        end
        
        
        function Pout = setdefault (Pin, varargin );
            % set the  default value for elements within a parameter set
            % Pin - input paramter array
            % varagin - name/value pairs to set
            Pout=Pin;
            % Take care of nested calls with variable number of parameters
            % where args appear as one element of a cell array
            
            while length(varargin)==1 & iscell (varargin{1}),
                varargin=varargin{1};
            end
            
            for i=1:2:length(varargin),
                j=pindex(Pin,varargin{i});
                Pout.paramDefault{j}=varargin{i+1};
            end
            
        end
        
        function Pout = setval (Pin, varargin );
            % set the value for a parameter set
            % Pin - input paramter array
            % varagin - name/value pairs to set
            Pout=Pin;
            % Take care of nested calls with variable number of parameters
            % where args appear as one element of a cell array
            
            while length(varargin)==1 & iscell (varargin{1}),
                varargin=varargin{1};
            end
            
            for i=1:2:length(varargin),
                name = varargin{i};
                value = varargin{i+1};
                j=pindex(Pin,name);
                if j ==0,
                    error (['Parameter:' name 'not found']);
                end
                if isempty(value),
                    disp([' Value not specified for Parameter:' name]);
                    disp('Parameter properties:');
                    dispFull(Pin(j));
                    Pout=Pin;
                    return
                end
                    
                if strcmp(Pin(j).paramType,'select'),
                    if ~any(strcmp(value,Pin(j).paramLimits)),
                        error (['The value: ''' value ''' is not in the limits list for select parameter: ' Pin(j).paramName ]);
                    end
                elseif strcmp(Pin(j).paramType,'number'),
                    if value > Pin(j).paramLimits(2),
                        error (['The value: ''' num2str(value) ''' too large for parameter: ' Pin(j).paramName ]);
                    elseif value < Pin(j).paramLimits(2),
                        
                    end
                 elseif strcmp(Pin(j).paramType,'logical'),
                     if ~islogical(value),
                          error (['Values for '  Pin(j).paramName ' must be type logical'  ]);
                     end
                end
                Pout(j).paramValue=value;
            end
        end
        
        
        function x = value(p);
            % Return value of a parameter
            if nargin==1,
                i=1;
            end
            if isnan(p.paramValue),
                x=p.paramDefault;
            else
                x=p.paramValue;
            end
            switch p.paramType
                case 'real'
                    x=double(x);
            end
            
        end
    end
end
% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt
