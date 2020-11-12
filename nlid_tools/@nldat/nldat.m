classdef nldat < nltop
    %nldat - data  class  for NLID toolbox
    % child of nltop.
    
    properties
        chanNames = {'x1' };
        chanUnits = { 'Default Units'  };
        domainIncr=1;
        domainName ='Time (s)';
        domainStart= 0;
        domainValues =NaN;
        dataSet = NaN;
        dataSize = NaN;
    end
    
    methods
        function d = nldat (z,varargin)
            if nargin ==0,
                % Generate an empty nldat variable
                return;
            elseif strcmp (class(z),'nldat');
                d =z;
                if nargin > 1,
                    set (d,[varargin(:)]);
                end
                
            elseif isa(z,'double')
                [nsamp,nchan,nreal]=size(z);
                if (nsamp == 1) & (nchan >10),
                    warning('Row vector. Transposing');
                end
                d.dataSet=z;
                [nsamp,nchan,nreal]=size(z);
                d.dataSize=size(z);
                for i=1:nchan,
                    d.chanNames{i} = [inputname(1) int2str(i) ];
                end
                
                if nargin > 1,
                    set (d,varargin{:});
                end
                
                
            elseif isa (z,'iddata');
                xin=z.InputData;
                xout=z.OutputData;
                d=cat(2,xin,xout);
                d=nldat(d);
                set (d,'domainIncr',z.Ts,'domainStart',z.Tstart);
            elseif isa (z,'parafcn');
                d=parafcn2nldat(z);
                if nargin > 1,
                    set (d,varargin);
                end
            elseif isa (z,'ms');
                M= full(z);
                mbyz=M.MbyZ{1};
                domainStart=mbyz(1);
                if length(mbyz)>1,
                    domainIncr=mbyz(2)-mbyz(1);
                else
                    domainIncr=nan;
                end
                d=nldat(M.Intensity{1},'domainStart',domainStart','domainIncr', ...
                    domainIncr, 'domainName','Da', 'comment',get(M,'comment'));
                
            elseif isa (z,'randv');
                d=randv2nldat(z);
                if nargin > 1,
                    set (d,varargin);
                end
            elseif isa (z,'signalv');
                d=signalv2nldat(z);
                if nargin > 1,
                    set (d,varargin);
                end
            else
                
                % Convert to a nldat object
                d=nldat;
                nlProps=properties(d);
                for i=1:length(nlProps),
                    set(d,nlProps{i},z.(nlProps{i}));
                end
            end
        end
        
        function sys = set.chanNames (sys, Value)
            if ~iscell(Value),
                sys.chanNames={Value};
            else
                sys.chanNames=Value;
            end
            
        end
        
        function sys=set.domainValues (sys, Value)
            
            if iscell(Value),
                sys.domainValues = Value;
            else
                sys.domainValues = double(Value);
            end
        end
        
        function sys = set.dataSet (sys, Value)
            sys.dataSet = Value;
            sys.dataSize=size(Value);
        end
        
        
        function sys = abs(sys)
            sys.dataSet=abs(sys.dataSet);
            sys.comment=[sys.comment  ';abs'];
            
        end
        
        function z = angle(x);
            z=x;
            z.dataSet=angle(x.dataSet);
            z.chanUnits='radians';
            z.comment=[x.comment  ';angle'];
        end
        
        function Z = cat (DIM, varargin);
            if nargin <3.
                error('cat must have at least three input arguments');
            end
            Z=(varargin{1});
            z=double(varargin{1});
            Z.comment = [ 'cat' int2str(DIM) ' ' inputname(2)  ];
            
            for i=2:nargin-1,
                y=nldat((varargin{i}));
                z=cat(DIM,z,double(y));
                Z.comment = [Z.comment ',' inputname(i+1)];
                switch DIM
                    case 2
                        Z.chanNames=cat(2,Z.chanNames,y.chanNames);
                end
            end
            Z.dataSet=z;
        end
        
        function [c,p] = corrcoef (x)
            % nldat wrapper for corrcoef function
            [nSamp,nChan,nReal]=size(x);
            xd=double(x);
            for i=1:nReal,
                zTemp=(xd(:,:,i));
                [cTemp,pTemp]=corrcoef(zTemp);
                if i==1,
                    c=cTemp;
                    p=pTemp;
                else
                    c=cat(3,c,cTemp);
                    p=cat(3,p,pTemp);
                end
            end
        end
        
        
        
        
        function z = chop (x, chopLen);
            % nldat/chop - chop data from begining and end of a record
            % z = chop (x, chopLen);
            % x - input data
            % chop len - length to chop off
            istart = round( (chopLen)/x.domainIncr) + 1 ;
            if istart>length(x.dataSet/2),
                error ('chopLen too large');
            end
            iend= length(x.dataSet) -istart + 1;
            x.dataSet=x.dataSet(istart:iend,:,:);
            x.domainStart=x.domainStart+chopLen;
            Comment = x.comment;
            set (x,'comment',[ Comment '; Chop']);
            z=x;
        end
        
        function z= cumtrapz (x)
            % nlid wrapper for cumtrapz
            % Z= cumtrapz (x)
            z=x;
            xd=double(x)*x.domainIncr;
            z.dataSet=cumtrapz(xd);
            z.comment=[ z.comment ' ; cumtrapz'];
        end
        
        
        function y = cumsum (x, DIM);
            
            y=nldat(x);
            if nargin==1,
                y.dataSet=cumsum((x.dataSet));
            else
                y.dataSet=cumsum(x.dataSet,DIM);
            end
            y.comment = [ 'CUMSUM of ' x.comment];
        end
        
        function [d_x] = delay(x,nDelay)
            %  del - delay the signal by n samples setting initial values to zero
            % x is the input signal
            % nDelay - dealy in samples
            d_x=x;
            [nSamp,nChan] = size(double(x));
            xStart=nldat(zeros(nDelay,nChan));
            xTemp=double(cat(1,xStart,x));
            xTemp=xTemp(1:nSamp,:);
            set(d_x,'dataSet',double(xTemp));
        end
        
        
        function v = ddt(x)
            % ddt function for nldat sets
            %	APPROXIMATES THE DERIVATIVE OF THE
            %	INPUT x BY THE FITTING A PARABOLA TO FIVE POINTS
            %	AND USING ITS ANALYTICAL DERIVATIVE.  THE EQUATION
            % 	FOR THIS APPROXIMATION IS:
            %
            %	    DX(N)=[-2X(N-2) - X(N-1) + X(N+1) +2X(N+2)]/(10*incr)
            %
            %	EXCEPT FOR THE FIRST TWO AND LAST TWO POINTS WHICH USE:
            %
            %	    DX(1)=   [-21X(1) + 13X(2) + 17X(3) - 9X(4)]/(20*INCR)
            %	    DX(2)=   [-11X(1) + 3X(2) + 7X(3) + X(4)]/(20*INCR)
            %	    DX(M-1)= [11X(M) - 3X(M-1) - 7X(M-2) - X(M-3)]/(20*INCR)
            %	    DX(M)=   [21X(M) - 13X(M-1) - 17X(M-2) + 9X(4M-3)]/(20*INCR)
            %
            
            [nsamp,nchan,nreal]=size(x);
            v=x;
            xd=x.dataSet;
            incr=x.domainIncr;
            for ireal=1:nreal,
                for ich=1:nchan,
                    vd(:,ich,ireal)=ddt_double(xd(:,ich,ireal),incr);
                end
            end
            v.dataSet=vd;
            v.comment =[ 'DDT ' x.comment ];
        end
        function [x] = del(x,delay)
            % this function adds a delay to the input signal
            % x is the input signal
            % delay - delay in (s)
            % d_x is the delayed singal
            % Now supports negative delays
            nDelay = floor(delay/get(x,'domainIncr'));
            xDataSet = get(x,'dataSet');
            d_x = zeros(size(x));
            if nDelay>=0
                d_x(nDelay+1:end,:) = xDataSet(1:end-nDelay,:);
            elseif nDelay<0
                nDelay = abs(nDelay);
                d_x(1:end-nDelay,:) = xDataSet(nDelay+1:end,:);
            end
            set(x,'dataSet',d_x);
        end
        function y = decimate (d, n, type);
            % usage  y = decimate (d, n);
            % or     y = decimate (d, n, 'fir');
            [nsamp,nchan,nreal]=size(d);
            y=nldat(d);
            if nargin == 3
                if (strcmp(type,'fir'))
                    for i=1:nchan,
                        for j=1:nreal,
                            xin = d.dataSet(:,i,j);
                            dout(:,i,j)=decimate(xin,n,'fir');
                        end
                    end
                else
                    error('Only argument is ''fir''')
                end
            else
                for i=1:nchan,
                    for j=1:nreal,
                        xin = d.dataSet(:,i,j);
                        dout(:,i,j)=decimate(xin,n);
                    end
                end
            end
            y.dataSet=dout;
            y.domainIncr=d.domainIncr*n;
            y.domainStart=d.domainStart + (n-1)*d.domainIncr;
            y.comment = ['DECIMATED ' d.comment ];
        end
        
        function Zout=detrend(ZIn)
            Zout=ZIn;
            [nsamp,nchan,nreal]=size(ZIn);
            for i=1:nchan,
                for j=1:nreal,
                    Zout.dataSet(:,i,j)=detrend(ZIn.dataSet(:,i,j));
                end
            end
            Zout.comment =[ ZIn.comment '; detrend'];
        end
        
        function d = domain(x);
            % returns domain vlaues for nldat object
            
            [nsamp,nchan,nreal]=size(x);
            if isnan(x.domainValues),
                d= ((1:nsamp)-1)*x.domainIncr + x.domainStart;
            else
                d=x.domainValues;
            end
            d=d(:);
        end
        
        
        function y=double(x)
            y = double (x.dataSet);
        end
        
        function i =end (A,k,n)
            i=builtin('end',A.dataSet,k,n);
        end
        
        function [P,pD]= fitdist(X,distributionName, binCenters)
            % Wrapper for matlab's fitdist  function
            % Usage: [P,pD]= fitdist(X,distributionName, binCenters)
            %   P - pdf object containing esimated and fited distributions.
            %   pD - matlabs ddistrobution object for the fit
            %   X - nldat object to compute distribttion (  1D for now);
            %   distributionName - same as for fitdsit
            % [binCenters] - point at which to compute distribution.
            %               Computed from range of X and nBins if not specified.
            
            x=X.dataSet;
            pD=fitdist(x,distributionName);
            if nargin <3
                P=pdf(X);
                binCenters=P.domainValues;
            else
                P=pdf(X,binCenters);
            end
            set(P,'chanNames','Estimate');
            y=pdf(pD,binCenters);
            P2=P;
            set(P2,'dataSet',y,'chanNames',[ distributionName ' fit'] );
            P=cat(2,P,P2);
            set(P,'comment', [ 'Estimated distribution and ' distributionName ' fit']);
            
        end
        
        
        
        
        
        function y = emean (x);
            % emean - ensemble meanif realizations
            
            [nsam, nchan, nresp]= size(x);
            y=nldat(x);
            xm=[];
            for i = 1:nchan
                xc=squeeze(x.dataSet(:,i,:))';
                xm= cat(2,xm, mean(xc)');
            end
            y.dataSet=xm;;
            set (y,'domainIncr',1,'comment','Ensemble Mean');
        end
        
        function z = ext (x, start, delt);
            % extract data based on domain values
            % x - input data
            % start - starting domain vlaue
            % delt  - length of extract. (default is to end of data set);
            % z = x(start:start_delt, :,:);
            
            istart = (start-x.domainStart)/x.domainIncr;
            if nargin == 2,
                iend=length(d.dataSet);
            else
                iend=min(istart + delt/x.domainIncr, length(x.dataSet));
            end
            x.dataSet=x.dataSet(istart:iend,:,:);
            x.domainStart=start;
            comment = get (x,'comment');
            set (x,'comment',[ comment '; ext']);
            z=x;
        end
        
        
        function y = extract (x, n1, n2);
            % Extract samples from a channel dropping first n1 and last n2 points
            if nargin==2,
                n2=n1;
            end
            [nsamp, nchan, nreal]=size (x);
            i1 = n1+1;
            i2 = nsamp - n2;
            y=x;
            dx=double(x);
            y.dataSet = dx(i1:i2,:,:);
            set (y,'comment','Extract');
        end
        
        function z = fft(x, nFFT);
            % fft function for nldat variables;
            % z = fft(x);
            % matlab's fft is applied to each realization of all channels in data set
            [nsamp, nchan, nreal]=size(x);
            if nargin ==1,
                nFFT=nsamp;
            end
            z=x;
            z.domainIncr= 1/(nFFT*x.domainIncr);
            set(z,'comment','FFT');
            z.domainStart=0;
            z.domainName='Frequency Hz';
            for ichan=1:nchan,
                for ireal=1:nreal,
                    z.dataSet(1:nFFT,ichan,ireal)=fft(x.dataSet(1:nFFT,ichan,ireal),nFFT);
                end
            end
        end
        
        function z = filter(x, B, A)
            % filter function for nldat variables;
            % z = filter(x, B,A);
            % filter is applied to each realization of all channels in data set
            if nargin <2
                error ('nldat/filt requires 2 or 3 inputs');
            end
            [nsamp, nchan, nreal]=size(x);
            z=x;
            set(z,'comment','filtered data ');
            for ichan=1:nchan,
                for ireal=1:nreal,
                    if nargin==2,
                        z.dataSet(1:end,ichan,ireal)=filter(B, x.dataSet(1:end,ichan,ireal));
                    elseif nargin==3,
                        z.dataSet(1:end,ichan,ireal)=filter(B,A, x.dataSet(1:end,ichan,ireal));
                    end
                end
            end
        end
        
        
        
        
        
        function Y = hwrect(X)
            % Half wave rectifier
            Y = X;
            x = double(X);
            y = (x>0).*x;
            Y.dataSet=y;
            Y.comment = [ X.comment '; hw rect'];
        end
        
        function xid = iddata (N)
            % converts nldata type to iddata type
            % only for single-input singple output
            
            xin=N.dataSet(:,1);
            xout=N.dataSet(:,2);
            xid=iddata(N.dataSet(:,2),N.dataSet(:,1));
            set(xid,'Ts',N.domainIncr');
        end
        
        function z = ifft(x, nFFT);
            % fft function for nldat variables;
            % z = ifft(x);
            % matlab's fft is applied to each realization of all channels in data set
            [nsamp, nchan, nreal]=size(x);
            if nargin ==1,
                nFFT=nsamp;
            end
            z=x;
            z.domainIncr= 1/(nFFT*x.domainIncr);
            set(z,'comment','IFFT');
            z.domainStart=0;
            z.domainName='Time s';
            for ichan=1:nchan,
                for ireal=1:nreal,
                    z.dataSet(1:nFFT,ichan,ireal)=ifft(x.dataSet(1:nFFT,ichan,ireal),nFFT);
                end
            end
        end
        
        
        function y = isnan (x);
            y=nldat(x);
            y=isnan(double(x));
        end
        
        function l = length (z)
            % length functon for nlad objects
            l = max(size(z));
        end
        
        function y = mean (x, DIM);
            if nargin==1,
                DIM=1;
            end
            y=nldat(mean(x.dataSet,DIM));
            set(y,'comment',['Mean(' inputname(1) ',' int2str(DIM) ')' ]);
            if DIM==3,
                set(y,'domainIncr',x.domainIncr,'domainName',x.domainName, 'domainStart',x.domainStart);
            end
        end
        
        function h=line(x,y)
            %  h=line(x,y,varargin)
            if nargin==1,
                % Nldat wrapper to line
                [nrow,ncol,nreal]=size(x);
                for i=1:ncol,
                    h(i)=line(domain(x),x.dataSet(:,ncol,1));
                end
            else
                h=line(double(x),double(y));
            end
        end
        
        
        
        
        
        function z = log(x)
            % Wrapper of log function for nldat objects
            z=x;
            set (z,'dataSet',log(x.dataSet));
        end
        
        function loglog(x,varargin)
            % wrapper for loglog
            arg=cat(2,{ 'xmode' 'log' 'ymode' 'log' } ,varargin{:});
            plot(x,arg);
        end
        
        
        function z = log10(x)
            % Wrapper of log10 function for nldat objects
            z=x;
            set (z,'dataSet',log10(x.dataSet));
        end
        
        function [z, iMax]  = max (x,y, DIM);
            x=nldat(x);
            [nSamp,nChan,nReal]=size(x);
            z=x;
            d=domain(x);
            if nargin==1,
                for iChan=1:nChan,
                    for iReal=1:nReal,
                        [xmax(iReal,iChan),iMax(iReal,iChan)]=max(x.dataSet(:,iChan,iReal));
                        domainValues(iReal,iChan,:)=d(iMax(iReal,iChan));
                    end
                end
                set(z,'dataSet',xmax,'domainValues',domainValues);
                
            elseif nargin==2,
                y=nldat(y);
                [ z.dataSet]=max(x.dataSet,y.dataSet);
            elseif nargin==3,
                if isempty(y)
                    [z.dataSet,iMax]=max(x.dataSet,[] ,DIM);
                else
                    y=nldat(y);
                    [z.dataSet]=max(x.dataSet,y.dataSet,DIM);
                end
            end
            set (z,'comment','Max');
        end
        
        function z = min (x,y, DIM);
            x=nldat(x);
            z=x;
            d=domain(x);
            [nSamp,nChan,nReal]=size(x);
            if nargin==1,
                for iChan=1:nChan,
                    for iReal=1:nReal,
                        [xmin(iReal,iChan),i]=min(x.dataSet(:,iChan,iReal));
                        domainValues(iReal,iChan,:)=d(i);
                    end
                end
                set(z,'dataSet',xmin,'domainValues',domainValues);
            elseif nargin==2,
                y=nldat(y);
                z.dataSet=min(x.dataSet,y.dataSet);
            elseif nargin==3,
                if isempty(y)
                    z.dataSet=min(x.dataSet,[] ,DIM);
                else
                    y=nldat(y);
                    z.dataSet=min(x.dataSet,y.dataSet,DIM);
                end
            end
            set (z,'comment','min');
        end
        
        function z = mod (x,y)
            z=x;
            set (z,'dataSet',mod(x.dataSet, double(y)));
        end
        
        function z = minus (x,y);
            % minus function for nldat variables;
            if isnumeric(y)
                y=nldat(y);
            end
            
            sx=size(x);
            xns =sx(1);
            xchan=sx(2);
            xreal=sx(3);
            sy= size(y);
            yns=sy(1);
            ychan=sy(2);
            yreal=sy(3);
            s1= [ 1 1 1];
            if sx==sy
                z=x;
                z.dataSet= x.dataSet - y.dataSet;
            elseif sx == s1
                z = y;
                z.dataSet= x.dataSet - y.dataSet;
            elseif sy == s1,
                z = x;
                z.dataSet= x.dataSet - y.dataSet;
            elseif  sx == [1 ychan 1 ]
                z=y;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) - y.dataSet(:,i,:);
                end
            elseif sy == [1 xchan 1 ];
                z=x;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) - y.dataSet(:,i,:);
                end
                % Separate value for each realization
            elseif sy == [1 xchan xreal ];
                z=x;
                for i=1:xchan,
                    for j=1:xreal,
                        z.dataSet(:,i,j) = x.dataSet(:,i,j) - y.dataSet(:,i,j);
                    end
                end
                % subtracting same value from all realizations
            elseif sy == [xns xchan 1 ];
                z=x;
                for i=1:xchan,
                    for j=1:xreal,
                        z.dataSet(:,i,j) = x.dataSet(:,i,j) - y.dataSet(:,i,1);
                    end
                end
            else
                error ('Dimension missmatch');
            end
            z.comment = [x.comment '-' y.comment];
        end
        
        function z = mtimes(x,y);
            % array multiply function for nldat variables;
            x=nldat(x);
            y=nldat(y);
            sx=size(x);
            xchan=sx(2);
            sy= size(y);
            ychan=sy(2);
            
            
            
            s1= [ 1 1 1];
            if sx==sy
                z=x;
                z.dataSet= x.dataSet .* y.dataSet;
            elseif sx == s1
                z = y;
                z.dataSet= x.dataSet .* y.dataSet;
            elseif sy == s1,
                z = x;
                z.dataSet= x.dataSet .* y.dataSet;
            elseif  sx == [1 ychan 1 ]
                z=y;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) .* y.dataSet(:,i,:);
                end
            elseif sy == [1 xchan 1 ];
                z=x;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) .* y.dataSet(:,i,:);
                end
            else
                error ('Dimension missmatch');
            end
            z.comment = [x.comment ' TIMES ' y.comment];
        end
        
        function z = phase(x);
            z=x;
            [nSap,nChan,nReal]=size(x);
            
            for i=1:nReal,
                f(:,1,i)=  unwrap(angle(x.dataSet(:,1,i)));
            end
            z.dataSet=f;
            z.comment=[ 'PHASE of ' x.comment];
        end
        
        function T = num2str(X,N)
            if nargin==1,
                T=num2str(double(X));
            else
                T=num2str(double(X),N);
            end
        end
        
        
        
        function z = lt (x,y);
            % lt function for nldat variables;
            if isnumeric(x)
                x=nldat(x);
            end
            if isnumeric(y)
                y=nldat(y);
            end
            
            sx=size(x);
            xchan=sx(2);
            sy= size(y);
            ychan=sy(2);
            s1= [ 1 1 1];
            if sx==sy
                z=x;
                z.dataSet= x.dataSet < y.dataSet;
            elseif sx == s1
                z = y;
                z.dataSet= x.dataSet < y.dataSet;
            elseif sy == s1,
                z = x;
                z.dataSet= x.dataSet < y.dataSet;
            elseif  sx == [1 ychan 1 ]
                z=y;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) < y.dataSet(:,i,:);
                end
            elseif sy == [1 xchan 1 ];
                z=x;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) < y.dataSet(:,i,:);
                end
            else
                error ('Dimension missmatch');
            end
            z.comment = [x.comment ' LT ' y.comment];
        end
        
        function z = plus (x,y);
            % plus function for nldat variables;
            if isnumeric(x)
                x=nldat(x);
            end
            if isnumeric(y)
                y=nldat(y);
            end
            
            sx=size(x);
            xchan=sx(2);
            sy= size(y);
            ychan=sy(2);
            s1= [ 1 1 1];
            if sx==sy
                z=x;
                z.dataSet= x.dataSet + y.dataSet;
            elseif sx == s1
                z = y;
                z.dataSet= x.dataSet + y.dataSet;
            elseif sy == s1,
                z = x;
                z.dataSet= x.dataSet + y.dataSet;
            elseif  sx == [1 ychan 1 ]
                z=y;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) + y.dataSet(:,i,:);
                end
            elseif sy == [1 xchan 1 ];
                z=x;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) + y.dataSet(:,i,:);
                end
            else
                error ('Dimension missmatch');
            end
            z.comment = [x.comment ' PLUS ' y.comment];
        end
        
        function z=power(x,n)
            z=x;
            z.dataSet=double(x).^n;
            z.comment=[z.comment ';.^' num2str(n)];
        end
        
        function s=prt (x)
            [nrow,ncol,nreal]= size(x);
            dv=domain(x);
            dt=get(x,'dataSet');
            cn=get(x,'chanNames');
            s=sprintf('%8s','-');
            for i=1:ncol,
                s=strcat( s,sprintf('%10s',char(cn{i})));
            end
            s1=[];
            for i=1:nreal,
                s1=sprintf ('%8s',char(dv(i)));
                for j= 1:ncol,
                    s1=strcat(s1,sprintf ('%10.4f', dt(1,j,i)));
                end
                s=strvcat(s,s1);
            end
        end
        
        function z=times(x,y)
            z=x;
            z.dataSet=x.dataSet.*y.dataSet;
            z.comment='x times y';
        end
        
        function z = real(x);
            z=x;
            z.dataSet=real(x.dataSet);
            z.comment = ['REAL part of ' x.comment];
        end
        
        function y = reshape (x, nsamp, nchan, nreal)
            % overlaid reshape for nldat objects
            if nargin < 4,
                nreal=1;
            end
            y=nldat(x);
            y.chanNames={};
            y.dataSet = reshape (x.dataSet, nsamp, nchan, nreal);
            y.dataSize=[ nsamp, nchan nreal];
            
            for i=1:nchan,
                y.chanNames{i}= ['Reshaped channel' num2str(i)];
            end
            y.comment = [' RESHAPED ' x.comment];
        end
        
        function y = reshapeChan (x, nSamp, nReal)
            %  reshape each channel in a nldat objects
            if nargin < 3,
                error ('3 reshapeChan requires three input parameters');
            end
            y=nldat(x);
            [nSampIn, nChan, nRealIn]=size(y);
            if nSamp*nReal ~= nSampIn*nRealIn
                error ('nsamp*nReal must be constant');
            end
            
            y.chanNames={};
            y.dataSet = [];
            y.dataSize=[ nSamp, nChan nReal];
            
            for iChan=1:nChan,
                yTemp=reshape (x.dataSet(:,iChan,:), nSamp,1,nReal);
                if iChan==1,
                    y.dataSet=yTemp;
                else
                    y.dataSet=cat(2, y.dataSet,yTemp);
                end
                y.chanNames{iChan}= ['Reshaped channel' num2str(iChan)];
            end
            y.comment = [' RESHAPED ' x.comment];
        end
        
        
        
        function z = rdivide(x,y);
            % array divide function for nldat variables;
            if isnumeric(y)
                y=nldat(y);
            end
            
            sx=size(x);
            xchan=sx(2);
            sy= size(y);
            ychan=sy(2);
            s1= [ 1 1 1];
            if sx==sy
                z=x;
                z.dataSet= x.dataSet ./ y.dataSet;
            elseif sx == s1
                z = y;
                z.dataSet= x.dataSet ./ y.dataSet;
            elseif sy == s1,
                z = x;
                z.dataSet= x.dataSet ./ y.dataSet;
            elseif  sx == [1 ychan 1 ]
                z=y;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) ./ y.dataSet(:,i,:);
                end
            elseif sy == [1 xchan 1 ] | sy == [ sx(1) sx(2) 1];
                z=x;
                for i=1:xchan,
                    z.dataSet(:,i,:) = x.dataSet(:,i,:) ./ y.dataSet(:,i,:);
                end
            elseif selse
                error ('Dimension missmatch');
            end
            z.comment = [ x.comment ' ./ ' y.comment];
        end
        
        function [Y, bestAC, bestVAF] = setAC ( X, ACin, nIterMax, vafLimit)
            % use stochastic interchange to generate a signal with the same
            % values but saoutocorrelation function that matchs AC in
            if nargin:4,
                vafLimit=99.5;
            end
            if nargin<3,
                nIterMax=1000; % number of unsuccesful interations
            end
            nIterTotalMax=10^6; 
                
            nSamp=length(X);
            nLags=get(ACin,'nLags');
            Y=X;
            bestAC=cor(Y, 'nLags',nLags); 
            bestVAF=double(vaf(ACin,bestAC,'total')); 
            jCount=1;
            for i=1:nIterTotalMax,
                iPerm=randperm(nSamp,2);
                newY=Y;
                newY.dataSet=Y.dataSet; 
                newY.dataSet(iPerm)=newY.dataSet(flip(iPerm));
                newAC=cor(newY, 'nLags',nLags);
                curVAF= double(vaf(ACin,newAC,'total'));
                jCount=jCount+1; 
                if curVAF>bestVAF,
                    %disp(['Best VAF: ' num2str(bestVAF) ' curVAF =' num2str(curVAF) ' jCount:' num2str(jCount)]);
                    Y=newY;
                    bestVAF=curVAF;
                    bestAC=newAC;
                    jCount=0;
                end
                if bestVAF>vafLimit,
                    disp('VAF limit reached'); 
                    break
                end
                if jCount>nIterMax,
                     disp('Unsuccessfull iteration limit reached'); 
                    break
                end          
            end
        end
        
        function [nsamp,nchan,nreal]= size (d, DIM)
            % overloaded size function for nldatclass.
            [nsamp,nchan,nreal]=size(d.dataSet);
            if nargout == 0 | nargout==1,
                if nargin ==1,
                    nsamp=[ nsamp nchan nreal];
                else
                    nsamp=size(d.dataSet, DIM);
                end
            elseif nargout ==2
                [nsamp,nchan]=size(d.dataSet);
            elseif nargout ==3
                return
            end
        end
        
        function y = smo(x,n, DIM);
            % smooth with a 3-point, zero-phase,  moving average filter
            %	y(i) = x(i-1)/4 + x(i)/2 + x(i+1)/4
            % USAGE	: y = smo(x,nSmooth)
            %	x	: input
            %	nSmooth: number of times x should be smoothed
            % See priovate function for detals
            
            if nargin==2,
                DIM=1;
            end
            y=x;
            xd=x.dataSet;
            yd=xd;
            [nsamp,nchan,nreal]=size(x);
            switch DIM
                case 1
                    
                    for ireal=1:nreal,
                        for ichan=1:nchan,
                            yd(:,ichan,ireal)=smo(xd(:,ichan,ireal),n);
                        end
                    end
                case 2
                    
                    for isamp=1:nsamp,
                        for ireal=1:nreal,
                            yd(isamp,:,ireal)=smo(xd(isamp,:,ireal),n);
                        end
                    end
                case 3
                    for isamp=1:nsamp,
                        for ichan=1:nchan,
                            yd(isamp,ichan,:)=smo(xd(isamp,ichan,:),n);
                        end
                    end
            end
            
            y.dataSet=yd;
            y.comment = [ 'SMOOTHED ' x.comment];
        end
        
        function y = squeeze(x);
            % squeeze for nldat objects
            y=x;
            D=squeeze(double(x));
            y.dataSet=D;
            y.comment=['SQUEEZE ' x.comment ];
        end
        
        function aOut=subsasgn(A,S,B)
            % subsasgn for nldat objects
            % V01-01 12 Nov 98
            
            
            if strcmp (S.type, '()')
                B=nldat(B);
                [nx(1),nx(2),nx(3)]=size(B.dataSet);
                for i=1:3,
                    ix{i}=1:nx(i);
                end
                for i=1:length(S.subs),
                    if ischar(S.subs{i}),
                        ix{i}=1:nx(i);
                    else
                        ix{i}=S.subs{i};
                    end
                end
                aOut=A;
                aOut.dataSet(ix{1},ix{2},ix{3})=B.dataSet;
                %
                % Make sure output channel names are defined
                [nSampOut, nChanOut, nRealOut]=size(aOut.dataSet);
                [nSampIn, nChanIn,nRealIn]=size(A);
                for iChan=1:nChanOut
                    if iChan <= nChanIn,
                        cnames(iChan)=A.chanNames(iChan) ;
                    else
                        cn1=cat(2,'x',int2str(iChan));
                        cnames=cat(2,cnames,{cn1});
                    end
                    aOut.chanNames=cnames;
                end
            elseif strcmp(S.type,'{}')
                error ('{} format not defined for nldat objects');
            elseif strcmp(S.type,'.')
                aOut=set(A,S.subs,B);
            end
        end
        
        
        
        function out = subsref (N, S)
            nTemp=N;
            for i=1:length(S),
                if strcmp(S(i).type,'.'),
                    nTemp=get(nTemp,S(i).subs);
                elseif strcmp(S(i).type,'()')
                    oDom=domain(nTemp);
                    d=nTemp.dataSet;
                    dTemp=builtin('subsref', d, S(i));
                    nTemp.dataSet=dTemp;
                    sIndex = S(i).subs;
                    % Fix DomainStart and DomainValues if necessary
                    if isnumeric(sIndex{1}),
                        nDom = oDom(sIndex{1});
                        diffIndex=diff(sIndex{1});
                        nTemp.domainStart=min(nDom);
                        if ~isnan(nTemp.domainValues) | length(unique(diffIndex))>1,
                            nTemp.domainValues=nDom;
                        elseif max(diffIndex)>1,
                            nTemp.domainIncr=nTemp.domainIncr*max(diffIndex);
                        end
                    end
                    % Fix Channel Names
                    if length(sIndex)>1,
                        nTemp.chanNames = nTemp.chanNames(S(i).subs{2});
                    end
                    %Fix for segdat
                    if isa(nTemp,'segdat')
                        onsetpointer = get(nTemp,'onsetPointer');
                        seglength = get(nTemp,'segLength');
                        set(nTemp,'onsetPointer', onsetpointer (S(i).subs{2}));
                        set(nTemp,'segLength' , seglength (S(i).subs{2}));
                    end
                end
                
                out = nTemp;
            end
        end
        
        function [xval,yval]= stairs (x)
            % options not yet implements
            [nsamp,nchan,nreal]=size(x);
            yval = x.dataSet;
            if isnan(x.domainValues),
                xval= ((1:nsamp)-1)*x.domainIncr + x.domainStart;
            else
                xval=x.domainValues;
            end
            if nargout==0,
                stairs (xval,yval);
            else
                [xval,yval]=stairs(xval,yval);
                %                     y=x;
                %                     set(y,'dataSet',y1,'domainValues',x1);
            end
        end
        
        
        function y = std (x, DIM);
            y=nldat(x);
            if nargin==1,
                DIM=1;
            end
            y.dataSet=std(x.dataSet,[],DIM);
            set(y,'comment',['std(' inputname(1) ',' int2str(DIM) ')' ]);
            
            if DIM==1,
                [nx,ny,nz]=size(y.dataSet);
                z=reshape(y.dataSet,ny,nz);
                y.dataSet=z';
                y.domainIncr=1;
                y.domainStart=1;
                y.domainName='Realization';
                for i=1:ny,
                    y.chanNames{i}=['std(' y.chanNames{i} ')'];
                end
            end
        end
        
        function y = sum (x, DIM);
            y=x;
            if nargin==1,
                y.dataSet=sum((x.dataSet));
            else
                y.dataSet=sum(x.dataSet,DIM);
            end
            y.comment=[ 'SUM of ' x.comment];
        end
        
        function z = uminus (x);
            % unary function for nldat variables;
            z=nldat(x);
            z.dataSet=-z.dataSet;
        end
        
        
        function z = uplus(x);
            % unary plus function for nldat variables;
            z=nldat(x);
        end
        
        function y = var (x);
            xd=double(x);
            [nSamp, nChan, nReal]=size(xd);
            for iChan=1:nChan,
                for iReal=1:nReal
                    y(iChan,iReal) =var(xd(:,iChan, iReal));
                end
            end
        end
        
        function errorPlot (d,e,  varargin)
            % nldat/plt - overloaded plot data with error values function for "nldat" class
            % plot (d, varagin);
            % options= 'plotmode' [line] 'plot type (line/xy/Super)
            %         'help_flag' 0 'display help (0=No/1=yes)
            %          nh - number of horizontal plots'
            %          nv - number of vertical plots'
            %          nplt - plot number'
            %          xmode
            %          y mode
            
            % Copyright 1999-2003, Robert E Kearney
            % This file is part of the nlid toolbox, and is released under the GNU
            % General Public License For details, see copying.txt and gpl.txt
            
            options={{'plotmode' 'line' 'plot type (line/xy/Super/bar)'} ...
                {'help_flag' 0 'display help (0=No/1=yes)'} ...
                {'line_color' '' 'Line color'} ...
                {'nh' NaN 'Number of horizontal plots'} ...
                {'nv' NaN 'Number of vertical plots'} ...
                {'nplt' 1 'Plot Number' } ...
                { 'xpanzoom' false 'enable linked zoom and pan for all axes'} ...
                { 'xpanwidth' NaN 'width of pan window' } ...
                {'xmode' 'linear' 'x axes mode [linear/log' } ...
                {'ymode' 'linear' 'y axes mode [linear/log/db' } ...
                
                };
            if isstr(d),
                arg_help('nldat/errorplot',options);
                return
            end
            if arg_parse_c('exact',options,varargin);
                return
            end
            %
            % Determine default layout
            %
            [nsamp, nchan, nreal]= size(d);
            if isnan (nv) & isnan (nh),
                % neither nv or nh set
                nh=ceil(nchan/2);
                nv=ceil(nchan/nh);
            elseif isnan (nv),
                % nh set
                nv=ceil(nchan/nh);
            elseif isnan(nh),
                % nv set
                nh = ceil (nchan/nv);
            end
            ;
            
            incr=d.domainIncr;
            start=d.domainStart;
            names=d.chanNames;
            if isnan(d.domainValues),
                t= ((1:nsamp)-1)*incr +start;
            else
                t=d.domainValues;
            end
            if strcmp(xmode,'log'),
                t=log10(t);
            end
            x=double(d);
            ed=double(e);
            if strcmp(ymode,'log'),
                x=log10(x);
                ed=log10(ed);
            end
            if strcmp(ymode,'db'),
                x=20*log10(x);
                ed=20*log10(d);
            end
            
            %
            % xy plot
            %
            % Line plots
            %
            if isnan(nplt),
                nplt=1;
            end
            
            for i=1:nchan,
                if (nchan > 1)
                    NP=nplt+i-1;
                    subplot (nv,nh,NP);
                end
                xd = squeeze(x(:,i,:));
                edCur=squeeze(ed(:,i,:));
                if nreal==1,
                    nlc = length(line_color);
                    if nlc>0,
                        plot (t,xd,line_color(max(nlc,i)));
                    else
                        plot (t,xd,'b', t, xd+edCur,'r', t, xd-edCur,'r');
                    end
                    
                    if strcmp(xmode,'log'),
                        xlabel([ 'Log ' d.domainName]);
                    else
                        xlabel(d.domainName);
                    end
                    switch ymode
                        case 'log'
                            ylabel([' log (' names{i} ')']);
                        case 'db'
                            ylabel([names{i} ' dB']);
                        otherwise
                            ylabel(names{i});
                    end
                end
                if xpanzoom,
                    xAxisPanZoom (gcf);
                end
                if ~isnan(xpanwidth),
                    xStart=get(d,'domainStart');
                    xEnd=xStart+xpanwidth;
                    set (gca,'xlim', [ xStart xEnd]);
                end
                title(d.comment);
            end
        end
    end
end

