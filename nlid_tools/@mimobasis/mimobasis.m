classdef mimobasis < nltop
    % mimobasis - multi-input, multi-output basis function expansion class for NLID toolbox.
    
    % mimobasis comprises of any 2-D basis function of input and scheduling variable (SV) 
    % whose coefficients are stored as polynomial objects in "coeffs" cell,
    % where each row i coeffs{i,:} represents an expansion with respect to SV for  
    % the i-th coefficient of the input expansion. For example, if SV is
    % expanded using Tchebychev polynomials, we will have: {tcheb0; tcheb1; … ; tcheb_n}
    % where n is the order of expansion with respect to input.
    
    % By Ehsan Sobhani, Release Notes: 
        % 1st Release on 2021-04-27 supports two kinds of 2D MIMO expansions: IRF and static nonlinerity (SNL)
    
    properties
        coeffs = {};
        nSVs = 1;
        nInputs = 1;
        svRange = NaN;
        inputRange = NaN;
        svExpType = '';        %-- Supported options are: {'tcheb'}
        svExpOrder = NaN;
        inputExpType = '';     %-- Supported options are: {'tcheb','laguerre'}
        inputExpOrder = NaN;
        coeffsStruct = NaN;
        parameterSet = param;
        notes = '2-D (input, SV) basis function expansion of two kind: IRF, Static NL';
    end
    
    methods
        %% Constructor of the class
        function basis = mimobasis(a, varargin)  
            if nargin == 0
                return
            elseif nargin == 1
                basis = nlmkobj(basis,a);
            elseif isa(a,'mimobasis')
                basis = nlmkobj(a,varargin{:});
            else
                basis = nlmkobj(basis,a,varargin{:});
            end
            %++ Type specific parameter: Decay parameter if input is expanded by Laguerre polynomials 
            switch basis.inputExpType
                case 'laguerre'
                    basis.parameterSet(1) = ...
                        param('paramName','alfa',...
                              'paramDefault',0.8,...
                              'paramHelp','Constant decay parameter of the Laguerre expansion',...
                              'paramType','number',...
                              'paramLimits',[0.01 1]); %-- The theoretical lower limit is any number bigger than 0
                otherwise 
                    return
            end           
        end
        
        %% % Overlaid set function for mimobasis objects
        function sys = set(sys, varargin)
            % v = varargin;
            if iscell(varargin{1}) && length(varargin{1})>1
                varargin = varargin{1};
            end
            
            if length(varargin)==1
                varargin = {varargin{1} []};
            end
            
            for i = 1:2:length(varargin)
                Prop = varargin{i};
                Value = varargin{i+1};
                % Check to see if Prop is the name of one of the properties of mimobase, and if yes ...
                % set it to the Value
                if ismember(Prop, properties(sys))
                    sys.(Prop) = Value;
                    if strcmp(Prop,'coeffs')
                        sys.inputExpOrder = length(Value)-1;
                    end
                    
                    if ~isempty(Value) && strcmp(Prop,'inputExpType')
                        switch Value
                            case 'tcheb'
                                outPs = param;
                            case 'laguerre'
                                % Laguerre needs a specific decay paramter (alfa) to define it
                                outPs(1) = param('paramName','alfa',...
                                           'paramDefault',0.8,...
                                           'paramHelp','Constant decay parameter of the Laguerre expansion',...
                                           'paramType','number',...
                                           'paramLimits',[0.01 1]); %-- The theoretical lower limit is any number bigger than 0
                            otherwise
                                outPs = param;    %-- only laguerre and tcheb are supported for inputExpType
                        end
                        sys.parameterSet = outPs;
                    end
                elseif ismember('parameterSet',fieldnames(sys)) %-- to check if 'parameterSet' is a property of mimobasis
                    % Check to see if the Prop is one of the parameterSets, and if yes ...
                    % the values of different parameters must change accordingly.
                    ps = sys.parameterSet;                      %-- Read the properties of the parameterSet object (from param class)
                    outPs = setval(ps,Prop,Value);
                    sys.parameterSet = outPs;
                end
            end
%             if ~isempty(inputname(1)),
%                 assignin('caller',inputname(1),sys);
%             end
        end
        
        %% Function to plot the MIMO Basis Function
        function plot(basis,varargin)
            options={{'n_bins_input' 40 'number of bins for a grid on input'} ...
                     {'n_bins_sv' 40 'number of bins for a grid on SV'} ...
            };
            if arg_parse(options,varargin)
                return
            end
            
            %++ Normalizing SV and generating SV grid
            if ~isempty(basis.svRange)
                %++ Calculate SV normalization parameters: Average and Range 
                maxSV = basis.svRange(2);
                minSV = basis.svRange(1);
                avg_sv = (maxSV + minSV) / 2;  %-- Note that this is not the statistical average
                rng_sv = maxSV - minSV;
                %++ Generate SV grid or operating points (OP)
                res_sv = rng_sv / n_bins_sv;
                sv_op = (minSV:res_sv:maxSV)';
                sv_op_n = (sv_op - avg_sv)*2 / rng_sv;
            else
                disp('The SV range is not specified in the object. Please set it and retry.')
                return
            end
            
            %++ Normalizing input and generating input grid
            if ~isempty(basis.inputRange)
                maxInput = basis.inputRange(2);
                minInput = basis.inputRange(1);
                rng_input = maxInput - minInput;
                res_input = rng_input / n_bins_input;
                input_op = (minInput:res_input:maxInput)';
                %++ Depending on the input expansion type, we may or may not need normalization
                switch basis.inputExpType
                    case 'laguerre'    %-- This case does not require normalization
                        input_op_n = input_op;
                    case 'tcheb'       %-- This case requires normalization
                        %++ Calculate input normalization parameter: Average  
                        avg_input = (maxInput + minInput) / 2;  %-- Note that this is not the statistical average
                        %++ Generate input grid or operating points (OP)
                        input_op_n = (input_op - avg_input)*2 / rng_input;
                    otherwise
                        disp('This input expansion type is not supported. Supported types are: laguerre and tcheb.')
                        return  
                end
            else
                disp('The input range is not specified in the object. Please set it and retry.')
                return
            end
            
            %++ Initialize mimo_surface and generate INPUT and SV grids
            nSV_op = length(sv_op_n);
            nInput_op = length(input_op_n);
            % mimo_surface = zeros(nSV_op,nInput_op);
            [INPUT,SV] = meshgrid(input_op_n,sv_op_n);
            
            %++ Converting the polynomial elements of the coeffs cell into a coefficient matrix named coeff_matrix
            q = basis.inputExpOrder;
            nRows = length(basis.coeffs);
            if nRows ~= q+1
                disp('The number of polynomials in the coeffs cell does not match the input expansion order!')
                return
            end
            
            %++ Future Extension: We assume all polynomials of SV have the same order, but ... 
            %++ in the future, that must change allowing different expansion orders w.r.t. SV
            sv_polynoms = basis.coeffs;                 %-- Reads all SV polynomial expansions
            m = sv_polynoms{1,1}.polyOrder;             %-- Read the order of the 1st SV polynomial {1,1}
            if m ~= basis.svExpOrder                    %-- In the Future Extension, we shall simply put a for-loop here for checking the order
                disp('The order of SV polynomials in the coeffs cell does not match the SV expansion order!')
                return
            end
            
            svPolyType = sv_polynoms{1,1}.polyType; 
            if ~strcmp(svPolyType,basis.svExpType)
                disp('The type of SV polynomials in the coeffs cell does not match the SV expansion type!')
                return
            end
                     
            %++ Generating the 2D basis function and ploting it as a 2D surface (or 3D plot)
            G_sv = multi_tcheb(sv_op_n,m);
            %++ Initialize and extract the coeff_matrix based on the verified orders 
            coeff_matrix = zeros(m+1,q+1);
            %++ This also depends on the input expnasion type
            switch basis.inputExpType
                case 'laguerre'
                    alfa = get(basis,'alfa');
                    G_input = laguerre(nInput_op,q,alfa);
                    for cnt_q = 1:q+1
                        coeff_matrix(:,cnt_q) = sv_polynoms{cnt_q,1}.polyCoef;
                    end
                
                case 'tcheb'
                    %++ Initialize and extract 
                    G_input = multi_tcheb(input_op_n,q);
                    for cnt_q = 1:q+1
                        coeff_matrix(:,cnt_q) = sv_polynoms{cnt_q,1}.polyCoef;
                    end
                
                otherwise
                    disp('This input expansion type is not supported. Supported types are: laguerre and tcheb.')
                    return
            end
            
            mimo_surface = G_sv * coeff_matrix * G_input';
            
            %++ Plotting the computed MIMO surface
            surf(INPUT,SV,mimo_surface,'LineStyle','none'); 
            xlabel('Input');  
            ylabel('SV'); 
            zlabel('Basis Amplitude');
            title('2D Basis Function');
        end
        
    end  %-- end of methods
end      %-- end of class definition

