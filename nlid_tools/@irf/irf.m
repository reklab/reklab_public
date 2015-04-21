classdef irf < kern
    % irf - impulse response function class for NLID toolboc.
    % parameters
    %   nLags - [< sample length] length of IRF in samples
    %   nSides -[1/2]  number of sides
    %  tvFlag- [Y/N] time varying IRF
    %  displayFlag -  [T/F] display plot
    %  irfErrorLvel -  [ real <=100] if confidence bounds are desired, this is the
    %                           computed
    %   irfIdMethod - method to estimate IRF
    %       tvfill - time varying identificiation
    %       cor - correlation method with Toeplitz
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
            I.parameterSet(j+1)=param('paramName','displayFlag','paramDefault','true',...
                'paramHelp','display');
            I.parameterSet(j+2)=param('paramName','irfErrorLevel','paramDefault',95, ...
                'paramHelp','error level','paramType','number', 'paramLimits',{0 100});
            I.parameterSet(j+3)=param('paramName','irfPseudoInvMode','paramDefault','full',...
                'paramHelp', 'pseudo-inverse order selection mode ', ...
                'paramType','select',...
                'paramLimits',{'full','auto','manual'});
            I.parameterSet(j+4)=param('paramName','irfIdMethod', ...
                'paramDefault','pseudo',...
                'paramHelp','Estimation method: tvfil/corr/pseudo/parmetric',...
                'paramType','select',...
                'paramLimits', {'tvfil', 'corr', 'param', 'pseudo'});
            
            I.comment='IRF Model';
            I.domainName='Lag (s)';
            
            if nargin==0;
                return
            elseif nargin==1,
                I=nlmkobj(I,a);
            elseif isa(a,'irf')
                I=nlmkobj(a,varargin{:})
            else
                I=nlmkobj(I,a,varargin{:});
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
            % CONSTRUCT an irf (impulse response function object
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
            [nsamp,nchan,nr]=size(z);
            
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
            % Parameter "tvFlag is "Yes" so Identify a time-varying impulse response functio
            %
            if strcmp(tvFlag,'Yes'),
                %  Paramaeter "estimationMethod" determines the method used:
                %   'tvfil' - finds least-squares solution using the data itself.
                %   'corr' - finds least-squares solution using correlation functions.
                %   'pseudo' - correlation function approach as with 'corr', but
                %   uses Dave's technique (adapted to tv case) to select singular vectors
                %   Default is 'corr'.
                x = squeeze(double(z(:,1,:)));
                y =squeeze(double(z(:,2,:)));
                conf_level=NaN;
                [Hident,bound,sing_vectors,cpu] = tv_ident(x,y,Ts,sides,numlags,irfIDMethod,conf_level);
                if nSides ==2,
                    LagStart=-nLags*Ts;
                    TimeStart=nLags*Ts;
                else
                    LagStart=0;
                    TimeStart=(nLags-1)*Ts;
                end
                
                set(i,'nSides',nSides, 'domainIncr',[Ts Ts],'dataSet',Hident, ...
                    'nLags',nLags,'domainStart',[TimeStart LagStart], 'domainName',{'Time' 'Lag'});
                i.kern=K;
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
                
                x = double(z(:,1,1));
                y =double(z(:,2,1));
                [filt,bounds]= fil_pinv(x,y,numlags,nSides, ...
                    irfErrorLevel,irfPseudoInvMode);
                if nr ==1,
                    dtotal=cat(2, filt/Ts, bounds/Ts);
                else
                    dreal=cat(2, filt/Ts, bounds/Ts);
                    dtotal=cat(3,dtotal,dreal);
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
    end
    
end

