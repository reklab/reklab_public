classdef fresp < nldat
    % frequency response class for NLID toolbox   
    
    properties
        parameterSet=param;
    end
    
    methods
        function F = fresp  (a,varargin)
            % Add object  specific parameters
            F.parameterSet(1)=param('paramName','Display_Flag','paramDefault','true',...
                'paramHelp','display','paramLimits', {'Yes' 'No'});
            F.parameterSet(2)=param('paramName','Method','paramDefault','tfe',...
                'paramHelp','Estimation method');
            F.parameterSet(3)=param('paramName','nFFT','paramDefault',NaN,...
                'paramHelp','segment length','paramType','number','paramLimits', [1 inf]);
            F.parameterSet(4)=param('paramName','nOverlap','paramDefault',0,...
                'paramHelp','Number of points window overlaps','paramType','number', 'paramLimits',[0 inf]);
            F.parameterSet(5)=param('paramName','Wind','paramDefault', NaN,...
                'paramHelp','Window ', 'paramType','number','paramLimits',[1 inf]);
            F.parameterSet(6)=param('paramName','Detrend_Mode','paramDefault','mean',...
                'paramHelp','Detrend (linear/mean/none' , ...
                'paramType','string', 'paramLimits',{ 'linear' 'mean' 'none'});
            
            I.comment='Frequency Response Model';
            if nargin==0;
                return
            elseif nargin==1,
                F=nlmkobj(F,a);
            elseif isa(a,'fresp')
                F=nlmkobj(a,varargin{:});
            else
                F=nlmkobj(F,a,varargin{:});
            end
            
        end
        
        function G = gain (F )
            %extract gain of fresp object as an nldat object
            G=nldat;
            G.dataSet=abs(F.dataSet(:,1,:));
            G.domainStart=F.domainStart;
            G.domainName=F.domainName;
            G.domainIncr=F.domainIncr;
            G.chanNames={'Gain'};
            G.comment='Frequency response gain';
        end
        
        function P = phase (F )
            %extract gain of fresp object as an nldat object
            [nSap,nChan,nReal]=size(F);
            P=phase(nldat(F)); 
            P.domainStart=F.domainStart;
            P.domainName=F.domainName;
            P.domainIncr=F.domainIncr;
            P.chanNames={'Phase'};
            P.comment='Frequency response phase';
        end
        
        function C = coherence (F )
            % extract gain of fresp object as an nldat object
            C=nldat;
            C.dataSet=abs(F.dataSet(:,2,:));
            C.domainStart=F.domainStart;
            C.domainName=F.domainName;
            C.domainIncr=F.domainIncr;
            C.chanNames={'coherence'};
            C.comment='Frequency response coherence';
        end
            
            
        
        function Fd = delay ( F, delt )
            % fresp/dealy  Add a delay to a fresp object
            
            d=2*pi*domain(F);
            delphi=d*delt;
            f=F.dataSet(:,1) ;
            gain = abs(f);
            phi = angle(f);
            newphi=phi - delphi;
            xreal=gain.*cos(newphi);
            ximag=gain.*sin(newphi);
            x=xreal + i*ximag;
            Fd=F;
            set(Fd,'dataSet',x);
        end
        
        
        function I = irf(F, nSides)
            if nargin==1,
                nSides=1;
            end
            I=irf;
            fd=F.dataSet;
            fpos=fd(:,1);
            lmax=length(fpos);
            fneg=conj(fpos(lmax:-1:2));
            f=[fpos;fneg];
            iout=real(ifft(f));
            Fincr= F.domainIncr;
            Tincr=1/(2*(Fincr*(lmax-1)));
            Tstart=0;
            nLags=length(iout);
            if nSides==2,
                iout=fftshift(iout);
                Tstart=-(length(iout)-1)*Tincr/2
                nLags=nLags/2;
            end
            
            set(I, 'dataSet',(iout/Tincr));
            set (I,'domainIncr',Tincr,'domainStart',Tstart);
            set (I, 'nLags', nLags);
        end
        
        function F = nlident (F, z,  varargin)
            % fresp/nlident
            % CONSTRUCT an fresp ( frequency response function object
            
            if nargin < 2,
                disp('NLIDtakes two inputs for fresp objects: irf, Z' );
            elseif nargin > 2,
                set(F,varargin);
            end
            % Parse options
            [nsamp,nchan,nr]=size(z);
            P=F.parameterSet;
            assign(P);
            ztype=class(z);
            % if spectrum
            switch ztype
                case 'irf'
                    error('obsolete now in irf'); 
                    F=irf2fresp(z);
                    fincr=z.domainIncr*length(z);
                    set (F,'domainName','Frequency','domainIncr',1/fincr);
                case 'nldat'
                    if isnan (nFFT)
                        nFFT=round(length(z)/10.);
                    end
                    if isnan (Wind),
                        Wind=hanning(nFFT);
                    end
                    Fs=1/z.domainIncr;
                    x=double(z(:,1));
                    y=double(z(:,2));
                    % [ t,f]=tfe (x,y,nFFT, Fs,Wind, NoOverlap,Detrend_Mode);
                    [ t,f]=tfestimate (x,y,Wind, nOverlap,nFFT, Fs);
                    [ c,f]=mscohere (x,y,Wind, nOverlap,nFFT, Fs);
                    %[ c,f]=cohere (x,y,nFFT, Fs, Wind, NoOverlap,Detrend_Mode);
                    set (F,'dataSet',cat(2,t,c));
                    set (F,'domainIncr', f(2)-f(1),'domainStart',f(1),'chanNames', {'tfe' 'cohere'}) ;
                    set(F,'domainName','frequency (Hz)');
                    set (F,'comment',['Frequency response:' z.comment]);
                otherwise
                    error (['fresp - does not support data type: ' ztype]);
            end
        end
        
        %
        %
        
        
        
        
        function y = nlsim ( F, xin )
            % Simulate response of fresp  to input data set
            % input options not fill defined as yet
            %
            
            % Copyright 2003, Robert E Kearney and David T Westwick
            % This file is part of the nlid toolbox, and is released under the GNU
            % General Public License For details, see ../copying.txt and ../gpl.txt
            
            xin=nldat(xin);
            f=F.dataSet;
            % strip out the coherence (second column of f)
            f = f(:,1);
            % generate negative frequency points in frequency response.
            % note that DC and fs/2 should only appear once, so remove the first
            % last points.
            fn=[f;conj(f(length(f)-1:-1:2))];
            B=ifft(fn);
            B=real(B(1:length(f)));
            %B = real(B);
            x=double (xin);
            yt=fftfilt(B,x,length(B));
            y=xin;
            set(y,'comment','filtered','dataSet',yt);
        end
        
        function plot (F,varargin)
            % plot a frequency response object
            %  xmode [linear/log] determines mode for plotting frequency axis
            options={{'mode' 'line' 'plot type (line/xy/Super/bar)'} ...
                {'help_flag' 0 'display help (0=No/1=yes)'} ...
                {'line_color' '' 'Line color'} ...
                {'nh' NaN 'Number of horizontal plots'} ...
                {'nv' NaN 'Number of vertical plots'} ...
                {'nplt' 1 'Plot Number' } ...
                {'moder' 'super' 'mode for plotting realizations' } ...
                {'offsetr' 1 'offset for plotting realizations' } ...
                {'xmode' 'linear' 'x axes mode [linear/log' } ...
                {'ymode' 'linear' 'y axes mode [linear/log/db' } ...
                };
            if isstr(F),
                arg_help('nldat/plt',options);
                return
            end
            if arg_parse(options,varargin);
                return
            end
            
             f=nldat(F); 
            f1=gain(F); %Gain
            f1=20*log10(double(f1)); 
            fPhase=phase(F); 
            f2=180*(unwrap(double(fPhase),[],1)/pi);
            [m,n,p]=size(F);
            if n==1,
                nv=2;
                f3=NaN;
            else
                nv=3;
                f3=f(:,2,:);
            end
            
            V=varargin;
            subplot (nv,1,1);
            f=nldat(F);
            title(' ');
            %
            % Get rid of 0's for log plots
            if strcmp(xmode,'log')
                imin=min(find(domain(f)>0));
                f=f(imin:length(f));
                f1=f1(imin:length(f1));
                f2=f2(imin:length(f2));
                f3=f3(imin:length(f3));
            end
            
            subplot (nv,1,1);
            set(f,'dataSet',f1,'chanNames',{'Gain'});
            plot(f,'xmode',xmode);
            ylabel('Gain (db)');xlabel('');
            subplot (nv,1,2);
            title(' ');
            set(f,'dataSet',f2,'chanNames',{'Phase'});
            plot(f,V);title('');
            ylabel('Phase (degrees)');
            if ~isnan(f3),
                xlabel('');
                subplot (nv,1,3);
                set(f,'dataSet',f3,'chanNames',{'Coherence'});
                plot(f,V);title('');
                set(gca,'ylim',[0 1.05]);
                title(' ');
            end
        end
        
        
        
        
    end
    
end

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt