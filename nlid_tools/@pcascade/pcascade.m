classdef pcascade < nlm
    % pcascade - parallel cascade model class for NLID toolbox.
    %  Parent: nlm
    properties
        nlmName='Parallel Cascade';
    end
    
    methods
        function PC = pcascade (z, varargin)
            % Initialize paramters
            PC.parameterSet(1)=param('paramName','idMethod', ...
                'paramDefault','eig', ...
                'paramHelp','Identification method' ,...
                'paramType','select', ...
                'paramLimits', {'eig' 'gen_eig' 'lm' 'slice' 'sls'});
            
            PC.parameterSet(2) = param('paramName','nPaths', ...
                'paramDefault',10,...
                'paramHelp' ,'Number of paths', ...
                'paramType','number',...
                'paramLimits',[1 1000]);
            PC.parameterSet(3)=param('paramName','nLags',...
                'paramDefault',16,...
                'paramHelp', 'Number of lags in each kernel' , ...
                'paramType','number', ...
                'paramLimits', [0 1000]);
            
            PC.parameterSet(4)=param('paramName','polyOrderMax',...
                'paramDefault',5, ...
                'paramHelp', 'Maximum order for series' , ...
                'paramType','number', ...
                'paramLimits', [0 10]);
            PC.parameterSet(5) = param('paramName','polyOrderSelectMode', ...
                'paramDefault','auto', ...
                'paramHelp','method to select polynomial order', ...
                'paramType','select',...
                'paramLimits',{ 'auto' 'full' 'manual' });
            
            PC.parameterSet(6) = param('paramName','nMaxPaths',...
                'paramDefault',10,...
                'paramHelp','Maximum number of paths tested in search',...
                'paramType','number', ...
                'paramLimits',[1 100]);
            PC.parameterSet(7) = param('paramName','testPaths', ...
                'paramDefault','yes',...
                'paramHelp','Retain only statistically significant paths',...
                'paramType','select', ...
                'paramLimits',{'yes','y','no','n'});
            PC.parameterSet(8) = param('paramName','nMaxIts', ...
                'paramDefault',100,...
                'paramHelp','Maximum number of paths tested in the expansion',...
                'paramType','number', ...
                'paramLimits',[1 10000]);
            PC.parameterSet(9)=param('paramName','searchThreshold', ...
                'paramDefault',.01, ...
                'paramHelp', 'NMSE for success', ...
                'paramType','number', ...
                'paramLimits',[0 1]);
            PC.parameterSet(10)=param('paramName','searchAccel', ...
                'paramDefault',.8, ...
                'paramHelp', 'ridge multiplied by accell after successful update',...
                'paramType','number', ...
                'paramLimits', [0.001 0.999]);
            PC.parameterSet(11)=param('paramName','searchDecel', ...
                'paramDefault',2, ...
                'paramHelp', 'ridge multipled by devel after unsuccessful update',...
                'paramType','number', ...
                'paramLimits', [1.0001 inf]);
            PC.parameterSet(12)=param('paramName','searchDelta', ...
                'paramDefault',10, ...
                'paramHelp', 'initial size of ridge added to Hessian',...
                'paramType','number',...
                'paramLimits',[0 inf]);
            PC.parameterSet(13)= param('paramName','initMethod', ...
                'paramDefault','random',...
                'paramHelp','Starting model specification', ...
                'paramType','select',...
                'paramLimits',{'current','random'});
            
            %%
            i=irf;
            t=polynom;
            set(t,'polyType','tcheb');
            elements = {i t; i t};
            set (PC,'elements',elements);
            set(PC,'idMethod','eig');
            if nargin==0;
                return
            elseif nargin==1,
                PC=nlmkobj(PC,z);
            else
                PC=nlmkobj(PC,z,varargin{:});
            end
        end
        
        
        
        function pc = nlident (pc,z,  varargin)
            % Estimate an pcascade - parallel cascade mode
            
            if nargin < 2,
                disp('nlident takes two inputs for pcascade objects: pcasacde, Z' );
            elseif nargin > 2,
                set(pc,varargin{:});
            end
            
            if isa(z,'nldat') | isa(z,'double')
                
                if isa(z,'nldat')
                    Ts=get(z,'domainIncr');
                else
                    subsys = get(pc,'elements');
                    f1 = pc{1,1};
                    Ts = get(f1,'domainIncr');
                    z = nldat(z,'domainIncr',Ts);
                end
                
               idMethod = lower(get(pc,'idMethod'));
                switch idMethod
                    case 'eig'
                        pc = pcm_eigen(pc,z,0);
                    case 'gen_eig'
                        pc = pcm_eigen(pc,z,1);
                    case 'slice'
                        pc = pcm_slice(pc,z);
                    case {'lm','sls'}
                        pc = pcm_nlls(pc,z);
                    otherwise
                        error('unsupported method');
                end
                % check mode, and proceed with identification
            elseif isa(z,'nlm')
                % convert z into a wiener-bose model
                ztype = lower(class(z));
                error('transformation from nlm to pcascade not yet implemented');
            else
                error('second argument must be a nldat object or a model with parent nlm');
            end
            set(pc,'comment',['Parallell cascade idMethod: ' idMethod]);
            
            
            
            function pc = pcm_eigen(pc,z,generalized)
                % parallel cascade identification via eigenvector method
                
                global PINV_ORDER;
                PINV_ORDER=0;
                
                u = z(:,1);
                y = z(:,2);
                yest = 0*y;
                zsize = size(z);
                N = zsize(1);
                hlen = get(pc,'nLags');
                
                P = getParamValStruct(pc.parameterSet);
                
                pcPath = lnbl;
                set(pcPath{1},'nLags',P.nLags);
                set(pcPath{2},'polyOrderMax',P.polyOrderMax, ...
                    'polyOrderSelectMode','auto');
                % try a first-order pathway
                
                set(pcPath,'initMethod','fil');
                resid = y - yest;
                ur = cat(2,u(hlen+1:N),resid(hlen+1:N));
                pcPath = nlident(pcPath,ur);
                  
                % check path for significance
                % if polynomial order = 0 minimizes the MDL, then
                % the linear dynamics are not significant.
                
                subsys = get(pcPath,'elements');
                h = subsys{1};
                m = subsys{2};
                
                order = get(m,'polyOrder');
                if order > 0
                    % pcPath is significant, so include
                    yest = nlsim(pcPath,u);
                    elements = {h m};
                    NPar = PINV_ORDER + order;
                    mdl_old = mdl_cost(y,yest,NPar,hlen);
                    disp(double(vaf(y,yest)));
                    num_paths = 1;
                else
                    elements = {};
                    NPar = 0;
                    mdl_old = mdl_cost(y,yest,NPar,hlen);
                    disp('no significant first-order path found');
                    num_paths = 0;
                end
                
                if generalized
                    set(pcPath,'initMethod','gen_eig');
                else
                    set(pcPath,'initMethod','eigen','idMethod','lm');
                   PINV_ORDER = hlen;
                end
                
                finished = 0;
                while ~finished
                    resid = y - yest;
                    ur = cat(2,u(hlen+1:N),resid(hlen+1:N));
                    pcPath = nlident(pcPath,ur);                 
                    subsys = get(pcPath,'elements');
                    h = subsys{1};
                    m = subsys{2};
                    order = get(m,'polyOrder');
                    if order < 1
                        finished = 1;
                    else
                        
                        pathout = nlsim(pcPath,u);
                        yest_new = yest + pathout;
                        NPar_new = NPar +PINV_ORDER + order -1;                        
                        mdl_new = mdl_cost(y,yest_new,NPar_new,hlen);
                        if mdl_new < mdl_old
                            elements = cat(1,elements,{h m});
                            yest = yest_new;
                            NPar = NPar_new;
                            mdl_old = mdl_new;
                            disp([ 'VAF:' num2str(double(vaf(y,yest))) ' MDL:' num2str(mdl_new)]);
                            num_paths = num_paths + 1;
                        else
                            finished = 1;
                        end
                    end
                    if num_paths == P.nMaxPaths;
                        finished=1;
                    end
                end
                
                set(pc,'elements',elements,'nPaths',num_paths);
            end
            
            function pc = pcm_slice(pc,z)
                % parallel cascade identification using single slices (original method).
                
                P = getParamValStruct(pc.parameterSet);
                Ts = get(z,'domainIncr');
                pcPath = lnbl;
                set(pcPath,'nLags',P.nLags,'polyOrderMax',P.polyOrderMax, ...
                    'idMethod','bussgang');
                subsys = get(pcPath,'elements');
                h = subsys{1};
                m = subsys{2};
                set(h,'domainIncr',Ts,'nLags',P.nLags);
                set(m,'polyOrderMax',P.polyOrderMax','polyType','tcheb');
                set(pcPath,'elements',{h m});
                u = z(:,1);
                y = z(:,2);
                resid = y;
                zsize = size(z);
                N = zsize(1);
                
                % Threshold for rejecting insingificant paths
                % expressed in tersm of Percent VAF -- hence 400 instead of 4
                % in the numerator
                if lower(P.testPaths(1)) == 'y'
                    TestSignificance = 1;
                    Threshold = 400/(N-P.nLags+1);
                else
                    TestSignificance = 0;
                end
                % Initialize the model
                num_paths = 0;
                elements = {};
                finished = 0;
                paths_tested = 0;
                
                while ~finished
                    path_order = ceil(3*rand(1));
                    switch path_order
                        case 1
                            set(pcPath,'lnInitMethod','correl');
                        case 2
                            set(pcPath,'lnInitMethod','slice2');
                        case 3
                            set(pcPath,'lnInitMethod','slice3');
                        otherwise
                            error('OOPS -- this error message should never get triggered');
                    end
                    ur = cat(2,u,resid);
                    ur0 = ur(P.nLags+1:N,:);
                    pcPath = nlident(pcPath,ur);
                    Significant = 1;
                    if TestSignificance
                        vf = double(vaf(pcPath,ur0));
                        if vf < Threshold
                            Significant = 0;
                        end
                    end
                    if Significant
                        % add the path to the model
                        num_paths = num_paths + 1;
                        elements = cat(1,elements,get(pcPath,'elements'));
                        set(pc,'elements',elements,'nPaths',num_paths);
                        yest = nlsim(pc,u);
                        resid = y - yest;
                    end
                    paths_tested = paths_tested+1;
                    if paths_tested > P.nMaxPaths
                        finished = 1;
                    end
                    if num_paths > P.nPaths
                        finished = 1;
                    end
                    disp(['nPaths: ' num2str(paths_tested) ' %VAF: ' num2str(double(vaf(pc,z)))]);
                end
            end
            
            function no_outputs = never_called(no_inputs)
                
                if nargin>2,
                    set(pc,varargin);
                end
                if isa(z,'lnl')
                    error('lnl2wiener conversion under development');
                end
                
                [nsamp,nchan,nreal]=size(z);
                numlags=get(pc,'nLags');
                if isnan(numlags),
                    numlags=min(64, nsamp/10);
                end
                zd=double(z);
                u = zd(:,1,:);
                y = zd(:,2,:);
                N = length(u);
                incr=get(z,'domainIncr');
                %
                % First Order path
                %
                resid = y;
                order=get(pc,'polyOrderMax');
                [h,m,outsum] = wiener_1(u,resid,numlags,order,'a');
                hirf=irf;
                set (hirf, 'dataSet',h,'nlags',numlags,'domainIncr',incr);
                mt= polynom;
                set (mt,'polyType','power');
                set (mt,'ployRange',m(1:2),'ployCoef',m(3:length(m)));
                elements = { hirf mt };
                if vaf(y,outsum) > 100*(numlags/N)
                    next_path = 2;
                    resid = y - outsum;
                    vaf_path = 100*(1 - (std(resid)/std(y))^2);
                    vafs(1,:) = [vaf_path vaf_path];
                    disp('Path 1 identified');
                    disp (['   VAF Path = ' num2str(vaf_path)]);
                else
                    next_path = 1;
                    disp('no significant first order path');
                    disp(vaf(y,outsum));
                end
                numpaths=get(pc,'nPaths');
                tol=get(pc,'PTolerance');
                for pathnum = next_path:numpaths
                    [h,m,out] = wiener_2(u,resid,numlags,order,'a');
                    hirf=irf; set (hirf, 'dataSet',h,'nLags',numlags,'domainIncr',incr);
                    mt= polynom;
                    set(mt,'polyType','tcheb');
                    set (mt,'polyRange',m(1:2),'polyCoef',m(3:length(m)));
                    vaf_path = 100*(1 - (std(resid-out)/std(resid))^2);
                    outsum = outsum + out;
                    resid = y - outsum;
                    vaf_total =  100*(1 - (std(resid)/std(y))^2);
                    vafs(pathnum,:) = [vaf_path vaf_total];
                    disp(['Path ' num2str(pathnum) ' identified']);
                    disp (['   VAF Path = ' num2str(vaf_path)]);
                    disp (['   VAF Total = ' num2str(vaf_total)]);
                    if vaf_path < tol,
                        disp (['VAF Path = ' num2str(vaf_path)]);
                        disp ([ 'VAF Path < tol. iteration terminated']);
                        pathnum=pathnum-1;
                        break
                    end
                    elements = cat (1, elements, { hirf mt });
                end
                set(pc,'clements',elements);
                set(pc,'nlags',numlags,'nPaths',pathnum);
            end
        end
    end
end

