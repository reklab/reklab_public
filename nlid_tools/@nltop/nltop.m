classdef nltop
    % nltop - parent class for all NLID classes
    
    properties
        comment
        version
    end
    
    methods
        function n = nltop(a, varargin)
            if nargin== 0,
                n.comment='Default comment';
                n.version='2.01';
            elseif isa (a,'nltop');
                n=a;
            else
                error(['nltop does not support inputs of class ' class(s)])
            end
        end
        
        function sys = set (sys, varargin)
            v=varargin;
            if iscell(varargin{1}) & length(varargin{1})>1,
                varargin=varargin{1};
            end
            if length(varargin)==1,
                varargin = { varargin{1} []};
            end
            for i=1:2:length(varargin),
                Prop=varargin{i};
                Value=varargin{i+1};
                % Check to see if it is a name is a property
                if ismember(Prop, properties(sys)),
                    sys.(Prop)=Value;
                elseif ismember('parameterSet',fieldnames(sys))
                    % Check to see if it is a parameter value
                    % Must change so that handles parameters of different
                    % names.
                    ps = sys.parameterSet;
                    outPs=setval(ps,Prop,Value);
                    sys.parameterSet=outPs;
                    %                     j=pindex(ps,Prop);
                    %                     if j,
                    %                         if isempty(Value), % Value not specified so show options
                    %                             dispFull(ps(j));
                    %                         else
                    %                             if strcmp(sys.parameterSet(j).paramType,'number'),
                    %                                 sys.parameterSet(j).paramValue={double(Value)};
                    %                             else
                    %                                 sys.parameterSet(j).paramValue={Value};
                    %                             end
                    %                         end
                    %                     else
                    %                         error([ Prop ' is not the name of a valid property or parameter']);
                    %                     end
                    %                 else
                    %                     error([ Prop ' is not the name of a valid property or parameter']);
                    %
                else
                    error([ Prop  ' is not a property or a paramrter']);
                end
            end
            if ~isempty(inputname(1)),
                assignin('caller',inputname(1),sys);
            end
        end
        
        function Value = get (sys, Property)
            % Get all public properties and their values
            if nargin==1
                Value= fieldnames(sys);
                if any(strcmp('parameterSet',Value)),
                    v=fieldnames(sys.parameterSet);
                    Value={ Value v};
                end
                
                
                % Handle various cases
            elseif nargin==2,
                % GET(SYS,'Property')
                if any(strcmp('parameterSet',fieldnames(sys))) & pindex(sys.parameterSet,Property),
                    Value=getParamValue(sys.parameterSet,Property);
                elseif any(strcmp(Property,fieldnames(sys))),
                    Value=sys.(Property);
                else
                    error (['Property not defined:' Property]);
                end
            end
            
            
        end
        
        function Value=subsref( sys, S )
            n=length(S);
            Value=sys;
            for iS=1:length(S),
                switch S(iS).type
                    case '.',
                        Value = get(Value, S(iS).subs);
                    case '{}'
                        nDim=length(S(iS).subs);
                        if nDim==1,
                            Value=sys.elements{S(iS).subs{1}};
                        elseif nDim==2,
                            Value=sys.elements{S(iS).subs{1},S(iS).subs{2}};
                        else
                            error ('Unexpected values');
                        end
                    otherwise
                        Value=builtin('subsref',Value,S(iS));
                end
            end
        end
        
        function Value= subsasgn (A,S,B)
            Value=A;
            sTemp=S;
            if length(sTemp)== 2,
                %disp('warning nlm/subsasgn may not work for this stucture');
                newA=subsref(A,sTemp(1));
                sTemp=sTemp(2);
                B= subsasgn( newA, sTemp, B);
            end
            switch S(1).type
                case '.'
                    Value = set(A,S(1).subs,B);
                case '{}'
                    i=S(1).subs{1};
                    if length(S(1).subs)==1,
                        j=1;
                    else
                    j=S(1).subs{2};
                    end
                    A.elements{i,j}=B;
                    Value=A;
                otherwise
                    Value = builtin('subsasgn',A,S,B);
                    
            end
        end
        
        
        
        
        
        function disp(sys)
            builtin('disp',sys);
            if ismember('parameterSet',properties(sys))
                disp('parameterSet:');
                disp(sys.parameterSet);
            end
            
        end
        
        function dispFull(sys)
            % Dispay object with full paramter information
            builtin('disp',sys);
            if ismember('parameterSet',properties(sys))
                disp('parameterSet:');
                dispFull(sys.parameterSet);
            end
            
        end
        
    end
end
% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt