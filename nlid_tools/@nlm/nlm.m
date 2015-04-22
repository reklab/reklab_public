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
            elseif isa(a,'nlm')
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
            subsys = sys.elements;
            [nparallel, nseries]=size(subsys);
            y=x(:,1)*0;
            for i=1:nparallel,
                xOut=x;
                for j=1:nseries,
                    ss=subsys{i,j};
                    xIn=xOut;
                    xOut = nlsim(ss,xIn);
                end
                y=y+xOut;
            end
            set(y,'comment','NLM simulation');
        end
        
        function plot (n)
            % Plot a nlm model
            e=n.elements;
            [nout,nin]=size(e);
            % SISO Series element
            [np,ns]=size(e);
            ifig=gcf;
            for i=1:np,
                for j=1:ns
                    subplot (np,ns,(i-1)*ns+j);
                    p=e{i,j};
                    plot (p);
                    %                     title(n.comment);
                    %                     h=get(gca,'title');
                    %                     u=get(h,'units');
                    %                     set(h,'units','pixels');
                    %                     p=get(h,'position');
                    %                     set(h,'position',p+[0 -30 0]);
                    %                     set(h,'units',u);
                    
                end
            end
            streamer ((n.comment));
        end
        
        
        
        
        function [npath,nel]= size (d)
            % overloaded size function for nlm.
            % npath - number of parallel paths
            % nel - number of elements per path
            E=d.elements;
            [npath, nel]=size(E);
            if nargout == 0,
                disp(['number of paths =' int2str(npath)]);
                disp(['number of series elements per path =' int2str(nel)]);
                
            elseif nargout == 1;
                npath =[ npath nel ];
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