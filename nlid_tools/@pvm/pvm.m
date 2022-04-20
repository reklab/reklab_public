classdef pvm < nltop
    % pvm - parameter-varying model parent class for NLID toolbox.
    
    % pvm model comprises series and/or parallel cascades of sub-models
    % stored in the elements property.
    % These submodels may be addressed using the format subModel= PVM{i,j}
    % where i is the path number and j the element in the path.
    % Properties of these submodel may be accessed/set using the form
    % PVM(i,j}.propertyName
    
    properties
        elements = { };
        inputName = { 'Input' };
        outputName = { 'Output' };
        schedVarName = { 'SV' };
        parameterSet=param;
        notes = 'PV model description / notes'
    end
    
    methods
        %% PVM Constructor: Constructs an instance of this class
        function PVM = pvm(a,varargin)
            if nargin==0
                return
            elseif nargin==1
                PVM = nlmkobj(PVM,a);
            elseif isa(a,'pvm')
                PVM = nlmkobj(a,varargin{:});
            else
                PVM = nlmkobj(PVM,a,varargin{:});
            end
        end
        
        %% PVM system elements display 
        function dispfull(sys)
            [m,n] = size(sys.elements);
            for i = 1:m
                for j = 1:n
                    disp(['element:' num2str(i) ',' num2str(j)]);
                    disp(sys.elements{i,j});
                end
            end
        end
        
        %% To Do: Add a method to get the pv system snapshots at a given array of SV values.
        
        %% PVM system simulator: Simulates the response of
        %++ various elements of pvm object to input and SV data
        %++ Each element has its own SV.
        function y = nlsim(sys,x,sv)
            subsys = sys.elements;
            [nparallel,nseries] = size(subsys);
            y = x(:,1)*0;
            for i = 1:nparallel
                xOut = x;
                for j = 1:nseries
                    ss = subsys{i,j};
                    xIn = xOut;
                    xOut = nlsim(ss,xIn,sv);
                end
                y = y + xOut;
            end
            set(y,'comment','PVM simulation');
        end
        
    end   %--> End of methods
end       %--> End of classdef

