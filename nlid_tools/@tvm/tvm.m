classdef tvm < nlm
    % Time Varyibg Model class
    % define and create time-varying model objects
    %  Parent: nlm
    
    % Copyright 2017, Robert E Kearney
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    
    properties
        tvStart=nan;
        tvIncr=0;
    end
    
    methods
        
        function TVM = tvm(a, varargin)
            TVM.comment='TV Model';
            TVM.parameterSet(1)=param('paramName','tvIdentMethod','paramDefault','manual',...
                'paramHelp', 'Method used to create tvm ', ...
                'paramType','select',...
                'paramLimits',{'manual','ensemble'});
            if nargin==0;
                return
            elseif nargin==1,
                TVM=nlmkobj(I,a);
            elseif isa(a,'nldat'),
                TVM = nlident(tvm, a, varargin{:});
            elseif isa(a,'irf')
                TVM=nlmkobj(a,varargin{:});
            else
                TVM=nlmkobj(I,a,varargin{:});
            end
        end
        
        function plot (TVM)
            el=TVM.elements;
            nEl=length(el);
            % Case for IRF models
            if isa(el{1},'irf'),
                
                for iEl=1:nEl,
                    curEl=el{iEl};
                    curVal=double(curEl);
                    if iEl==1,
                        Z=zeros(nEl,length(curVal));
                    end
                    Z(iEl,:)=curVal;
                    
                end
                X=domain(TVM);
                Y=domain(curEl);
                mesh(Y,X,Z);
                xlabel('Lag');
                ylabel('Time');
                title(TVM.comment);
                % Case for polynomial models
            elseif isa(el{1},'polynom')
                for iEl=1:nEl,
                    curP=el{iEl};
                    if iEl==1,
                        pRange=curP.polyRange;
                        xVal=linspace(pRange(1),pRange(2));
                        Z=zeros(nEl,length(xVal)) ;
                    end
                    Z(iEl,:)=double(nlsim(curP,xVal'));
                end
                
                X=domain(TVM);
                mesh(xVal,X,Z);
                xlabel('Input Value');
                ylabel('Time');
                zlabel('Output Value');
                title(TVM.comment);
            elseif isa(el{1},'nlm')
                firstEl=el{1};
                [nRow,nCol]=size(firstEl);
                figNum=0;
                nFig=nRow*nCol;
                for iRow=1:nRow,
                    for iCol=1:nCol
                        curBlock=block(TVM,iRow,iCol);
                        figNum=figNum+1;
                        subplot (nRow,nCol,figNum);
                        modelType=class(curBlock.elements{1});
                        set(curBlock,'comment', [ curBlock.comment ': ' modelType]);
                        plot(curBlock);
                    end
                end
                
            else
                error('model type not supported');
                
            end
            
        end
        
        function TVM = nlident ( TVM, Z, varargin)
            disp(Z);
            modelPrototype=varargin{1};
            modelType =class(varargin{1});
            switch modelType
                case 'irf'
                    disp('estimate a TV irf');
                    X=squeeze(double(Z(:,1,:)));
                    Y=squeeze(double(Z(:,2,:)));
                    dt=Z.domainIncr;
                    nSides=modelPrototype.nSides;
                    nLags=modelPrototype.nLags;
                    approach=modelPrototype.irfIdMethod;
                    conf_level=modelPrototype.irfErrorLevel;
                    if strcmp(approach,'pseudo'),
                        [hIdent,bound,sing_vectors,cpu] = tvIdentEnsemble(X,Y,dt,nSides,nLags,approach,conf_level);
                    else
                        [hIdent] = tvIdentEnsemble(X,Y,dt,nSides,nLags,approach);
                    end
                    if nSides ==2,
                        LagStart=-nLags*dt;
                        TimeStart=nLags*dt;
                    else
                        LagStart=0;
                        TimeStart=(nLags-1)*dt;
                    end
                    i=modelPrototype;
                    set(i,'nSides',nSides, 'domainIncr',dt, ...
                        'nLags',nLags,'domainStart', LagStart, 'domainName','Lag');
                    
                    [nReal,nLag]=size(hIdent);
                    I={};
                    for iReal=1:nReal,
                        curIRF=hIdent(iReal,:);
                        set(i,'dataSet', curIRF);
                        I{iReal}=i;
                    end
                    set(TVM,'tvStart',TimeStart,'tvIncr',dt,'elements',I');
                    
                otherwise
                    disp(['tvm: identification not supported for model  type: ' modelType ]);
            end
            %% Reformat for output.
            
            
            disp(modelType);
        end
        
        
        
        function yPre = nlsim( TVM, xIn)
            % Simulate the ensemble response of a TVM model
            if TVM.tvIncr ~= xIn.domainIncr
                error ('TV Model and input must have the same sampling rate');
            end
            indexStart=1 +(TVM.tvStart-xIn.domainStart)/TVM.tvIncr;
            index=indexStart;
            dt=TVM.tvIncr;
            yPre=xIn;
            [m,n,p]=size(xIn);
            yData=zeros(m,n,p);
            el=TVM.elements;
            nEl= length(TVM.elements);
            for iEl=1:nEl,
                curModel=el{iEl};
                if isa(curModel,'irf'),
                    yTemp=  tvIRFsim(curModel.dataSet,xIn.dataSet,index,curModel.nSides, dt);
                    yData(index,1,:)=yTemp;
                elseif isa(curModel,'polynom')
                    yTemp=nlsim(curModel,xIn(index,1,:));
                    yData(index,1,:)=double(yTemp);
                elseif isa(curModel,'nlbl')
                    pCur=curModel{1};
                    xTemp=xIn;
                    xTemp=nlsim(pCur,xIn);
                    curIRF=curModel{2};
                    yTemp= tvIRFsim(curIRF.dataSet,xTemp.dataSet,index,curIRF.nSides,dt );
                    yData(index,1,:)=yTemp;
                    
                else
                    error ('TV Model not supported');l
                end
                index=index+1;
            end
            set(yPre,'dataSet',yData);
            
            
            
            function y = tvIRFsim(H, X, i, nSides, dt)
                % Simulate the response of a TV IRF
                [n,n,r]=size(X);
                hlen = length(H);
                z=zeros(hlen,1);
                i=i+hlen;
                for iReal=1:r
                    Xpad=cat(1,z,X(:,1,iReal),z);
                    
                    if nSides == 2
                        x = Xpad(i-(hlen-1)/2:i+(hlen-1)/2,:);
                    else
                        x = Xpad(i-hlen+1:i,:);
                    end
                    y(iReal) = (flipud(H(:))'*x*dt)';
                end
            end
            
            function y = tvNLBLsim(NLBL, X, i)
                % Simulate the response of a TV NLBL
                p=NLBL{1};
                hIRF=NLBL{2};
                nSides=hIRF.nSides;
                H=double(hIRF); H=H(:);
                dt=X.domainIncr;
                hlen = length(H);
                if nSides == 2
                    x = X(i-(hlen-1)/2:i+(hlen-1)/2,1,:);
                else
                    x = X(i-hlen+1:i,1,:);
                end
                xp=nlsim(p,x);
                x=double(xp); x=squeeze(x);
                y = (flipud(H)'*x*dt)';
            end
        end
        
        
        function NN= normGainLE ( N )
            % Normalize a TV nlblmodel  to assign all gain to the NL element.
            % Scale Hammerstein model apropriately
            % Only valid for low pass systems.
            NN=N;
            nElements=length(N.elements);
            for iElement=1:nElements,
                curElement=N.elements{iElement};
                if ~strcmp(class(curElement),'nlbl'),
                    error('Model must be of type nlbl');
                end
                normElement=normGainLE(curElement);
                NN.elements{iElement}=normElement;
                set(NN,'comment','Normalized nlbl');
            end
        end
        
        function tvBlock = block (TVM, row, col );
            % Return block row,col from TVM
            tvBlock=TVM;
            el=TVM.elements;
            nEl=length(TVM.elements);
            subElements={};
            for iEl=1:nEl,
                curEl=el{iEl};
                tempEl=curEl{row,col};
                subElements{iEl}=tempEl;
            end
            tvBlock.elements=subElements;
            
        end
        
        function d = domain(TVM)
            %Return domain values for a tvm
            nEl=length(TVM.elements);
            dt=TVM.tvIncr;
            t=TVM.tvStart-dt;
            for i=1:nEl
                t=t+dt;
                d(i)=t;
            end
        end
        
        function Value=subsref( sys, S )
            %% {} returns a cell array with individual models for  specified indices
            %% () returns a TVM with elements specified, as  ell array of more thn model
            n=length(S);
            Value=sys;
            for iS=1:length(S),
                switch S(iS).type
                    case '.',
                        Value = get(Value, S(iS).subs);
                    case '{}'
                        nDim=length(S(iS).subs);
                        if nDim==1,
                            Value={Value.elements{S(iS).subs{1}}};
                            if length(Value)==1,
                                Value=Value{1};
                            end
                        elseif nDim==2,
                            Value=Value.elements{S(iS).subs{1},S(iS).subs{2}};
                        else
                            error ('Unexpected values');
                        end
                    case '()'
                        nDim=length(S(iS).subs);
                        if nDim==1,
                            t=domain(Value);
                            t=t(S(iS).subs{1});
                            Value.tvStart=t(1);
                            Value.elements={Value.elements{S(iS).subs{1}}};
                        elseif nDim==2,
                            Value.elements=Value.elements{S(iS).subs{1},S(iS).subs{2}};
                        else
                            error ('Unexpected values');
                        end
                    otherwise
                        Value=builtin('subsref',Value,S(iS));
                end
            end
        end
        
        
    end
end

