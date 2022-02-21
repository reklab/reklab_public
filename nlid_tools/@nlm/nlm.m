classdef nlm < nltop
    % nlm - nonlinear model parent class for NLID toolbox.
    
    % nlm model comprise series and or parallel cascades  of sub-models
    % stored in the elements property.
    % These submodels may be addressed using the format subModel= NLM{i,j}
    % where i s th epath number and j the element in the path.
    % Properties of these submodel  may be access/set using the form
    % NLM(i,j}.propertyName
    
    
    properties
        elements = { };
        inputName = { 'Input' };
        outputName = { 'Output' };
        parameterSet=param;
        notes = 'Notes'
        
    end
    
    methods
        function N = nlm (a,varargin)
            % Add object  specific parameters
            N.comment='IRF Model';
            if nargin==0;
                return
            elseif nargin==1,
                N=nlmkobj(N,a);
            elseif isa(a,'nl')
                N=nlmkobj(a,varargin{:});
            else
                N.nlmkobj(N,a,varargin{:});
            end
        end
        
        function dispfull(sys)
            [m,n]=size(sys.elements);
            for i=1:m,
                for j=1:n
                    disp(['element:' num2str(i) ',' num2str(j)]);
                    disp(sys.elements{i,j});
                end
            end
        end
        
        
        
        function y = double (x);
            % Double for nlm objects
            y= cat(1, x.range, x.coef);
        end
        
        function disp(sys)
            builtin('disp',sys);
            if ismember('parameterSet',properties(sys))
                disp('parameterSet:');
                disp(sys.parameterSet);
            end
            
        end
        
        function y = nlsim ( sys, x )
            % Simulate response of nlm object to input data set
            % Handle segdat objects separately to simulate properly. 
            if isa(x,'segdat')
                [nSamp,nChan,nReal]=size(x);
                if nChan ~=2,
                    error ('nlsim needs ID output data for segdat simulations');
                end
                xOutId=x(:,2);
                
            end
            subsys = sys.elements;
            [nparallel, nseries]=size(subsys);
            y=x(:,1)*0;
            for i=1:nparallel,
                xOut=x;
                for j=1:nseries,
                    ss=subsys{i,j};
                    if isa(x,'segdat') & isa(ss,'ssm')
                        xIn=cat(2,xOut,xOutId);
                    else
                        xIn=xOut;
                    end
                    xOut = nlsim(ss,xIn);
                end
                y=y+xOut;
            end
            set(y,'comment','NLM simulation');
        end
        
        function plot (n, nV, nH,nSub)
            %  plot an nlm moodel: plot (n, nv, nh, nSub)
            % n - nm model
            % optional parameters (must specify all or nonee
            % nv - number of vertical panels in plot
            % nh - number of horiontal pannels in plot 
            % nSub - subplot por elements in n
              e=n.elements;
            [nout,nin]=size(e);
            nElement=nout+nin-1;
            switch nargin
                case 1
                    nV=nout;
                    nH=nin;
                    nSub=1:nElement;
                case 4
                    if length(nSub) ~=nElement
                        error ('nSub muust equal the number of elements in the model');
                    end
                otherwise
                   error ('nlm.plot requires either 1 ot 4 inputs');
            end
          
            
            % SISO Series element
            [np,ns]=size(e);
            ifig=gcf;
            for i=1:np,
                for j=1:ns
                    subplot (nV,nH,nSub((i-1)*ns+j));
                    p=e{i,j};
                    p.comment= [p.comment ' of ' n.comment];
                    plot (p);
                   
                end
            end
            
        end
        
        
        
        
        function [npath,nel, nReal]= size (d)
            % overloaded size function for nlm.
            % npath - number of parallel paths
            % nel - number of elements per path
            % nReal is faked from first element 
            E=d.elements;
            [m,n,nReal]=size(d.elements{1,1});
            [npath, nel, nReal]=size(E);
            if nargout == 0,
                disp(['number of paths =' int2str(npath)]);
                disp(['number of series elements per path =' int2str(nel)]);
                
            elseif nargout == 1;
                npath =[ npath nel nReal];
            end
        end
        
       
        
        
        function vf = vaf (model,data,spare)
            % VAF function for NLM class models
            %
            % vf = vaf(model,data)
            %
            % where model is any model of class nlm
            % and data is an nldat object containing [input,output]
            
            if nargin == 1
                error('I/O data must be supplied if first argument is of class NLM');
            end
            
            if nargin == 2
                % data contains both input and output
                if isa(data,'nldat');
                    ddata = double(data);
                    Ts = data.domainIncr;
                elseif isa(data,'double')
                    ddata = data;
                    Ts = 0;
                else
                    error('second argument must be either double or nldat');
                end
                [N,m] = size(ddata);
                if min(N,m)==1
                    error('must supply both input and output to compute vaf');
                elseif m > 2
                    error('To test mimo systems, vaf(model,inputs,outputs)');
                end
                u = nldat(ddata(:,1));
                y = nldat(ddata(:,2));
                if Ts > 0
                    set(u,'domainIncr',Ts);
                    set(y,'domainIncr',Ts);
                end
            end
            
            if nargin == 3
                if isa(data,'nldat');
                    u = data;
                elseif isa(data,'double');
                    u = nldat(data);
                    Ts = 0;
                else
                    error('second argument must be either double or nldat');
                end
                if isa(spare,'nldat');
                    y = spare;
                    if Ts == 0
                        Ts = get(y,'domainIncr');
                        set(u,'domainIncr',Ts);
                    end
                elseif isa(data,'double');
                    y = nldat(spare);
                    if Ts > 0
                        set(y,'domainIncr',Ts);
                    end
                else
                    error('third argument must be either double or nldat');
                end
            end
            
            yest = nlsim(model,u);
            if Ts == 0
                Ts = get(yest,'domainIncr');
                set(y,'domainIncr',Ts);
            end
            vf = vaf(y,yest);
        end
    end
    
end