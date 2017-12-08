classdef tvm < nlm
    % Time Varyibg Model class
    %  Parent: nlm
    % Time-varying model is stored as a cell array of nlm models with one
    % element for each time.
    %
    % nlsim method supports nlmparallel cascade  models of arbitrary
    % complexity comprising parallel cascades of IRF and POLYNOM elements.
    %
    % nlident method supports identification of:
    %   irf, polynom and nl block models.
    %
    % tvIdentMethod is a parameter defining how the tv model was created
    % Values are:
    %   manual - model created manually
    %   ensemble - model estimated using an ensemble idnetification method
    %   basisexpansion - model estimate using a basis expansion mehtod
    %
    % see tvmDemo for demonstration of use.
    %
    %
    %Copyright 2017, Robert E Kearney
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    % need to docment how to use.
    properties
        tvStart=0;
        tvIncr=1;
    end
    
    methods
        
        function TVM = tvm(a, varargin)
            TVM.comment='TV Model';
            TVM.parameterSet(1)=param('paramName','tvIdentMethod','paramDefault','manual',...
                'paramHelp', 'Method used to create tvm ', ...
                'paramType','select',...
                'paramLimits',{'manual','ensemble' 'basisexpansion'});
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
        
        function tvmIdent = nlident ( TVM, Z, modelPrototype, varargin)
            % tvIdent = nlident ( tvm, Data, modelProtype,
            % Identify a tvm model
            % Inuts:
            % TVM - template for the tv model
            %   TVM.tcvIdentMethod determines identification method
            %       = ensemble - ensemble TV idnetifiction method
            %       = basisexpansion = temporal expansion method
            %  Z  - input/output data set
            %  modelPrototype - protype of TV model ( irf, polynom)
            %   - propeties of the prototype determine those of the model
            %   estimates.
            % varargin - name/value pais defining options for temporal
            % expansion identifiction.
            % vrargin(1) = BASIS - polynom object defining the basis function used
            %         withi the temporal expansion method
            %        'periodic' 'no' 'Periodic data {yes/no}'
            %         'method' 'Bayes' 'Linear estimation algorithm to be used {Bayes/OLS}'}...
            modelType =class(modelPrototype);
            switch modelType
                case 'irf'
                    tvmIdent = nlidentIRF ( TVM, Z, modelPrototype,  varargin{:});
                case 'polynom'
                    tvmIdent = nlidentPOLYNOM ( TVM, Z, modelPrototype );
                case 'nlbl' % Hammerstein system 
                     tvmIdent = nlidentNLBL ( TVM, Z, modelPrototype,  varargin{:});
                otherwise
                    disp(['tvm: identification not supported for model  type: ' modelType ]);
            end
            
        end
        
        function tvmIdent = nlidentIRF ( TVM, Z, modelPrototype, varargin )
            % Identify a TV IRF Model
            if nargin>3,
                BASIS=varargin{1};
                options={{'periodic' 'no' 'Periodic data {yes/no}'}...
                    {'method' 'Bayes' 'Linear estimation algorithm to be used {Bayes/OLS}'}...
                    };
                if arg_parse_c('exact',options,varargin(2:end));
                    return
                end
            end
            tvmIdent=TVM;
            X=squeeze(double(Z(:,1,:)));
            Y=squeeze(double(Z(:,2,:)));
            dt=Z.domainIncr;
            t=domain(Z);
            nSides=modelPrototype.nSides;
            nLags=modelPrototype.nLags;
            approach=modelPrototype.irfIdMethod;
            conf_level=modelPrototype.irfErrorLevel;
            assign(TVM.parameterSet)
            switch tvIdentMethod
                case 'ensemble',
                    
                    if strcmp(approach,'pseudo'),
                        [hIdent,bound,sing_vectors,cpu] = tvIdentEnsemble(X,Y,dt,nSides,nLags,approach,conf_level);
                    else
                        [hIdent] = tvIdentEnsemble(X,Y,dt,nSides,nLags,approach);
                    end
                    if nSides ==2,
                        LagStart=-nLags*dt;
                        TimeStart=Z.domainStart+nLags*dt;
                    else
                        LagStart=0;
                        TimeStart=Z.domainStart+(nLags-1)*dt;
                    end
                case 'basisexpansion'
                    B=basisfunction(BASIS,t); B=double(B);
                    [hIdent, x_pred, Extra] = tv_irf_ident_expansion(X, Y, B, nLags, nSides, dt,periodic, method);
                    if nSides ==2,
                        LagStart=-nLags*dt;
                    else
                        LagStart=0;
                    end
                    TimeStart=Z.domainStart;
                otherwise
                    error 'Bad Value for tvIdentMethod');
            end
            
            i=modelPrototype;
            set(i,'nSides',nSides, 'domainIncr',dt, ...
                'nLags',nLags,'domainStart', LagStart, 'domainName','Lag');
            
            [nReal,nLag]=size(hIdent);
            I={};
            for iReal=1:nReal,
                curIRF=hIdent(iReal,:);
                set(i,'dataSet', curIRF(:));
                I{iReal}=i;
            end
            set(tvmIdent,'tvStart',TimeStart,'tvIncr',dt,'elements',I');
        end
        
        
        
        function tvmIdent = nlidentPOLYNOM ( TVM, Z, modelPrototype )
            % identify a TV polynomial from an ensemble
            tvmIdent=TVM;
            [nSamp, nChan, nReal]=size(Z);
            P=[];
            for iSamp=1:nSamp
                z1=squeeze(Z(iSamp,1,:));
                z2=squeeze(Z(iSamp,2,:));
                z=cat(2,z1,z2);
                pTemp=nlident(modelPrototype,z);
                P{iSamp}=pTemp;
            end
            set(tvmIdent,'tvStart',Z.domainStart,'tvIncr',Z.domainIncr,'elements',P');
        end
        
        function tvmIdent = nlidentNLBL ( TVM, Z, modelPrototype, varargin )
            % Identify a TVNLBL Model
            
            tvmIdent=TVM;
            X=squeeze(double(Z(:,1,:)));
            Y=squeeze(double(Z(:,2,:)));
            dt=Z.domainIncr;
            polyPrototype=modelPrototype{1};
            irfPrototype=modelPrototype{2};
            polyOrder=polyPrototype.polyOrderMax
            nSides=irfPrototype.nSides;
            nLags=irfPrototype.nLags;
            if nSides==1
                M1=0;
                M2=nLags-1;
            else
                 M1=round(nLags/nSides);
                 M2=M1;
            end          
            tol = .01;
               [coeffs,range,IRF,num_iterations] = tv_nlbl_ident_ensemble(X,Y,dt,polyOrder,M1,M2,tol);

           
            assign(TVM.parameterSet)
%             switch tvIdentMethod
%                 case 'ensemble'
%                     [coeffs,range,IRF,num_iterations] = tv_nlbl_ident_ensemble(X,Y,dt,polyOrder,M1,M2,tol);
                    
                    if nSides ==2
                        LagStart=-nLags*dt;
                        TimeStart=Z.domainStart+nLags*dt;
                    else
                        LagStart=0;
                        TimeStart=Z.domainStart+(nLags-1)*dt;
                    end

            
            i=irfPrototype;
            set(i,'nSides',nSides, 'domainIncr',dt, ...
                'nLags',nLags,'domainStart', LagStart, 'domainName','Lag');
            p=polyPrototype;
            n=modelPrototype;
            [nReal,nLag]=size(IRF);
            I={};P={};
            N={}; 
            for iReal=1:nReal,
                curIRF=IRF(iReal,:);
                set(i,'dataSet', curIRF(:));
                I{iReal}=i;
                set(p,'polyCoef',coeffs{iReal});
                P{iReal}=p;              
                set(n,'elements',{p i});
                N{iReal}=n;
            end
            
            set(tvmIdent,'tvStart',TimeStart,'tvIncr',dt,'elements',N');
        end
        
        
        function V=tvResid ( tvModel, Z);
            simY=nlsim(tvModel,Z(:,1,:));
            Y=Z(:,2,:);
            subplot(3,1,1);
            plot(Y);
            title('Output');
            subplot (3,1,2);
            plot (simY);
            title('Estimated Output');
            subplot(3,1,3);
            plot(Y-simY);
            V=vaf(Y,simY,'total');
            title(['Residuals VAF= ' num2str(V)]);
            
        end
        
        
        function yPre = nlsim( TVM, xIn)
            % yPre= nlsim( TVM, xIn) simulate a TVM model
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
                % Simulate the response of a TV for one element of an ensemble
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
        
        function nlmtst(sys)
            tvmDemo;
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

