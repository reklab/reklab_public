classdef irf < kern
    % irf - impulse response function class for NLID toolboc.
    % parameters
    %   nLags - [< sample length] length of IRF in samples
    %   nSides -[1/2]  number of sides
    %  tvFlag- [true/false] time varying IRF
    %    TV IRFs are stored using realizations for the time variable: (lag,1,time);
    %  tvStartTime - stsart time for TV ITfs
    %  displayFlag -  [T/F] display plot
    %  irfErrorLvel -  [ real <=100] if confidence bounds are desired, this is the
    %                           computed
    %  irfIdMethod - method to estimate IRF
    %       tvfill - time varying identificiation
    %       corr - correlation method with Toeplitz
    %       pseudo - pseduo inverse
    %  irPseudoInvMode - [ full/auto/manual] mode for pseudo inverse method
    %       full - retain all singular values
    %       auto - determine how many to retain automatically based on MDL
    %       manual - select umber of singular values auotomatically
    % Copyright 2003, Robert E Kearney and David T Westwick
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    
    properties
        Name='irf';
        irfBounds=NaN;
    end
    
    methods
        function I = irf (a,varargin)
            % Add IRF specific parameters
            j=length(I.parameterSet);
            I.parameterSet(j+1)=param('paramName','tvFlag','paramDefault',false, ...
                'paramType','logical', ...
                'paramHelp','tvFlag( true false) kernel is time varying');
            I.parameterSet(j+2)=param('paramName','tvStartTime','paramDefault',0, ...
                'paramType','number', ...
                'paramHelp','start time for TV parameter');
            I.parameterSet(j+3)=param('paramName','displayFlag','paramDefault','true',...
                'paramHelp','display');
            I.parameterSet(j+4)=param('paramName','irfErrorLevel','paramDefault',95, ...
                'paramHelp','error level','paramType','number', 'paramLimits',[0 100]);
            I.parameterSet(j+5)=param('paramName','irfPseudoInvMode','paramDefault','auto',...
                'paramHelp', 'pseudo-inverse order selection mode ', ...
                'paramType','select',...
                'paramLimits',{'full','auto','manual'});
            I.parameterSet(j+6)=param('paramName','irfIdMethod', ...
                'paramDefault','pseudo',...
                'paramHelp','Estimation method: tvfil/corr/pseudo/param/user',...
                'paramType','select',...
                'paramLimits', {'tvfil', 'corr', 'param', 'pseudo' 'user'});
            I.parameterSet(j+7)=param('paramName','irfFigNum', ...
                'paramDefault',1,...
                'paramHelp','Figure to use for manual selection [1]',...
                'paramType','number',...
                'paramLimits', [1 inf]);
            
            I.comment='IRF Model';
            I.domainName='Lag (s)';
            
            if nargin==0;
                return
            elseif nargin==1,
                I=nlmkobj(I,a);
            elseif isa(a,'irf')
                I=nlmkobj(a,varargin{:});
            else
                I=nlmkobj(I,a,varargin{:});
            end
            
        end
        
        function t = domain(I),
            nLags=length(I.dataSet);
            nSides=get(I,'nSides');
            nLags=get(I,'nLags');
            if nSides==1,
                N=nLags;
                startT=0;
            else
                N=nLags*2+1;
               startT=-nLags*I.domainIncr; 
            end
            for i=1:N,
                t(i)=startT+ (i-1)*I.domainIncr;
            end
            
        end
        
        
        function F = fresp(I)
            % Convert an IRF to frequency response
            F=fresp;
            k=I.dataSet;
            fr=fft(k);
            l=1+(floor(length(fr)/2));
            fr=fr(1:l);
            Tincr=I.domainIncr;
            fftincr=1/(Tincr*(length(k)-1));
            fr=fr*Tincr;
            set(F,'dataSet',fr,'domainStart',0, 'domainIncr',fftincr, ...
                'domainName','Hz');
        end
        
        function [beta, irfEst]=fitKBI (I, beta0)
            % I -= irf to fit
            % beta0 - initial parameter estimates [ k b i];
            % beta - estimated parameters and properties
            % irfEst - estimated IRF function
            x=domain(I)';
            y=double(I);
            beta=fitnlm(x,y,@fkbiFnc, beta0)
            irfTmp=fkbiFnc(beta.Coefficients.Estimate,x);
            irfEst=I;
            set(irfEst,'dataSet',irfTmp,'comment','IRF estimate');   
            clf; 
            plot(I);
            h=line(irfEst);
            set(h,'color','r');
        end

        function I = irf2 (I, varargin  )
            % Generate a second order IRF with gain, damping and natural
            % freuquency
            paramList={{'g' 1 'gain'} ...
                {'z' .75 'Damping '} ...
                {'w' 12 'Natural Frequency'} ...
                {'delt' .01 'Time increment '} ...
                {'irflen', 1 'IRF length'} ...
                };
            
            if arg_parse(paramList,varargin);
                return
            end
            t=0:delt:irflen;
            i=fgzw(t, [ g z w]);
            set (I,'domainIncr',delt,'nLags',length(t),'dataSet', i, ...
                'comment','Second order irf');
        end
        
        function l = length(I)
            nSides=get(I,'nSides');
            nLags=get(I,'nLags');
            
            if nSides==1
                l=nLags;
            else
                l=nLags*2 + 1;
            end
        end
        function I=nlmtst(i)
            irfDemo
        end
        
        
        function plotBounds(I)
            %- plot IRF with error bounds
            plot(I);
            hold on
            t=domain(I);
            v=I.dataSet;
            b=I.irfBounds;
            plot(t,v+b,'r',t,v-b,'r');
            hold off
            C=I.comment;
            assign(I.parameterSet);
            title([C ' with ' num2str(errorLevel) '% error bounds'])
            
        end
        
        function i = nlident (i, z,  varargin)
            % identify an irf (impulse response function object
            % Setup default values
            %
            
            if nargin < 2,
                disp('NLIDtakes two inputs for irf objects: irf, Z' );
            elseif nargin > 2,
                set(i,varargin{:});
            end
            if isa(z,'fresp'),
                P=get(i,'parametersSet');
                assign(P);
                i=fresp2irf(z, nSides);
                return
            end
            if isa(z,'double'),
                [nrow,ncol]=size(z);
                if (nrow <2),
                    z=z';
                    [nrow,ncol]=size(z);
                    warning('irf/nlident:Input has too few rows. Transposing and continuing');
                end
                if ncol==1,
                    %
                    % one Channel so it is the IRF
                    K=i.kern;
                    set(K,'dataSet',z);
                    i.kern=K;
                    return
                elseif ncol==2,
                    %
                    % Two Channels so do identification
                    %
                    z=nldat(z);
                else
                    %
                    %MOre than 2 channel is an error
                    error('Input data cannot have more than 2 channels');
                end
            end
            
            
            
            % Parse options
            [nsamp,nchan,nReal]=size(z);
            %% Nchan=1 so convert it to an IRF
            if nchan==1,
                i.chanNames=z.chanNames;
                i.chanUnits=z.chanUnits;
                i.domainIncr=z.domainIncr;
                i.domainName=z.domainName;
                i.domainStart=z.domainStart;
                i.domainValues=z.domainValues;
                i.dataSet=z.dataSet;
                i.dataSize=z.dataSize;
                i.comment=z.comment;
                return
            end
            
            
            assign(i.parameterSet);
            if isnan (nLags)
                nLags=min(64, round(nsamp/100));
            end
            Ts=z.domainIncr;
            if isnan(Ts),
                error('IRF estimation requires domainIncr to be specified');
            end
            
            %
            % force two sided IRFs to be odd length
            %
            if nSides==2,
                numlags=2*nLags+1;
                sides='two';
            else
                numlags=nLags;
                sides='one';
            end
            %
            % Parameter "tvFlag is true so Identify a time-varying impulse
            % response from ensemble data
            %
            if (tvFlag),
                %  Paramaeter "estimationMethod" determines the method used:
                %   'tvfil' - finds least-squares solution using the data itself.
                %   'corr' - finds least-squares solution using correlation functions.
                %   'pseudo' - correlation function approach as with 'corr', but
                %   uses Dave's technique (adapted to tv case) to select singular vectors
                %   Default is 'corr'.
                x = squeeze(double(z(:,1,:)));
                y =squeeze(double(z(:,2,:)));
                [m,n]=size(x);
                zeroPad=zeros(nLags, n);
                x=cat(1,zeroPad, x);
                y=cat(1,zeroPad,y);
                if nSides==2,
                    x=cat(1,x,zeroPad);
                    y=cat(1,y,zeroPad);
                end
                
                
                conf_level=NaN;
                [hIdent,bound,sing_vectors,cpu] = tv_ident(x,y,Ts,sides,numlags,irfIdMethod,conf_level);
                if nSides ==2,
                    LagStart=-nLags*Ts;
                    TimeStart=nLags*Ts;
                else
                    LagStart=0;
                    TimeStart=(nLags-1)*Ts;
                end
                [m,n]=size(hIdent);
                tvI=hIdent';
                tvI=reshape(tvI,n,1,m);
                set(i,'nSides',nSides, 'domainIncr',Ts,'dataSet',tvI, ...
                    'nLags',nLags,'domainStart', LagStart, 'domainName','Lag', 'tvStartTime',z.domainStart);
                %
                % Otherwise identify a time-invariant impulse response function
                %
            else
                %    Parmeter "mode" determines the method used to determin the rank of the
                %   identifiable subspace used in the low-rank projection
                %   'full' - full length
                %   'auto' - automatically determine
                %   'manual' - determine length interactively
                %               if 0<mode<numlags+1, then mode is used as the pseudo-inverse
                %               order.
                for iReal=1:nReal,
                    zd=double(z);
                    x = double(zd(:,1,iReal));
                    y =double(zd(:,2,iReal));
                    [filt,bounds]= fil_pinv(x,y,numlags,nSides, ...
                        irfErrorLevel,irfPseudoInvMode, irfFigNum);
                    if iReal ==1,
                        dtotal=cat(2, filt/Ts, bounds/Ts);
                    else
                        dreal=cat(2, filt/Ts, bounds/Ts);
                        dtotal=cat(3,dtotal,dreal);
                    end
                end
                
                if nSides==1,
                    start=0;
                elseif nSides==2,
                    start=-(numlags-1)*Ts/2;
                end
                
                set(i,'nSides',nSides, ...
                    'domainIncr',Ts,'dataSet',dtotal(:,1,:), 'nLags',nLags,'domainStart',start, ...
                    'domainName',z.domainName);
                i.irfBounds=dtotal(:,2,:);
            end
            
        end
        
        function y = nlsim ( model, xin )
            % irf/nlsim Simulate response of IRF to input data set
            % input options not fully defined as yet
            filter = model.dataSet;
            if isa(xin,'double'),
                xin=nldat(xin);
                set(xin,'domainIncr',model.domainIncr);
            end
            delx = xin.domainIncr;
            deli=model.domainIncr;
            if delx ~= deli,
                W=(str2mat('Model & data have different domain increments', ...
                    'the output of the IRF depends on the sampling rate', ...
                    'Output may be scaled incorrectly and/or have the wrong increment'));
                warning(' ');disp(W)
            end
            x=double (xin);
            
            incr = model.domainIncr;
            assign(model.parameterSet);
            %
            % Simulate a time-varying response
            %
            [irfLen, irfDim, nSampIrf]=size(model);
            [nSamp, nChan, nReal]=size(xin);
            
            if (tvFlag),
                if nSamp>nSampIrf,
                    warning(' Number of input samples > number of TV IRfs. Output truncated' );
                    nSamp=nSampIrf;
                end
                if xin.domainStart ~= tvStartTime,
                    error('Data and TV IRFs have different start times');
                end
                
                if nSides==1,
                    offSet=1;
                else
                    offSet= -(model.domainStart/model.domainIncr) + 1;
                end
                irf=double(model);
                for iSamp=1:nSamp,
                    yTemp=0;
                    for jLag=1:irfLen,
                        k=iSamp -jLag +offSet;
                        if (k>=1 & k<=nSamp),
                            yTemp=yTemp+(irf(jLag,1,iSamp)*x(k,1,:));
                        end
                    end
                    yOut(iSamp,1,:)=yTemp*incr;
                end
                y=xin;
                set(y,'comment','filtered','dataSet',yOut);
                %
                % Simulate a time-invariant response
                %
            else
               % x=x(:,1,:);
                [nsamp,nchan,nreal]=size(x);
                [nsampFilter, nchanFilter,nrealFilter]= size(filter);
                % check that realizations are reasonable
                if nreal==nrealFilter | nreal==1 | nrealFilter==1,
                    nRealMax=max(nreal,nrealFilter);
                else
                    error ('Number of realizatins in IRF and input are not compatible');
                end
                for i=1:nchan,
                    for j=1:nRealMax,
                        jFilter=min(j,nrealFilter);
                        jInput=min(j,nreal);
                        yout(:,i,j) = filter_ts(filter(:,jFilter), x(:,i,jInput), nSides, incr);
                    end
                end
                y=xin;
                set(y,'comment','filtered','dataSet',yout);
                
            end
            
        end
        
        
    end
end


