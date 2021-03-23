classdef spect < nldat
    % spect - power spectrum class for NLID toolbox. 
    
    properties
        
       parameterSet=param;
    end
    
    methods
        function s = spect  (a,varargin)
            % Add object  specific parameters
            
            s.parameterSet(1) =param('paramName','nFFT','paramDefault',NaN, ...
                'paramHelp','Length of FFT',  ...
                'paramType','number',...
                'paramLimits',[1 inf]);
            s.parameterSet(2)=param('paramName','nOverlap', ...
                'paramDefault',NaN,...
                'paramType','number',...
                'paramHelp','Overlap between FFT segments');
            s.parameterSet(3)=param('paramName' ,'windowFunction', ...
                'paramDefault',NaN, ...
                'paramHelp','Window function', ...
                 'paramType','number'...
                 );
            s.parameterSet(4)=param('paramName','spectType', ...
                'paramDefault','Auto', ...
                'paramHelp','Type of spectrum', ...
                'paramType','string'  ...
                );
            
            s.comment='Spectrum of ';
            if nargin==0;
                return
            elseif nargin==1,
                s=nlmkobj(s,a);
            elseif isa(a,'spect')
                s=nlmkobj(a,varargin{:});
            else
                s=nlmkobj(s,a,varargin{:});
            end
            
        end
        %-------------------------------------------------------------------------------------
        function S = nlident (Sin, z, varargin)
            % identify spect objects
            
            S=Sin;
            z=nldat(z);
            if nargin < 2,
                disp('nlident takes two inputs for cor objects: cor, Z' );
            elseif nargin > 2,
                set(S,varargin);
            end
            
            if any(isnan(double(z))),
                error('Signal contains NaNs');
            end
            tempComment=S.comment;
            
            [nsamp,nchan,nreal]=size(z);
            incr=z.domainIncr;
            ninputs=1;
            noutputs=nchan-ninputs;
            
            assign(S.parameterSet);
            if isnan(nFFT),
                nFFT=nsamp/10;
            end
            nFFT = round(nFFT);
            if isnan(nOverlap)
               nOverlap=0;
            end
            if isnan(windowFunction),
                windFunction=nFFT;
            end
            windowFunction=windFunction';
            incr=z.domainIncr;
            fs=1/incr;
            Order=1;
            % noutputs =0 so we have an autospectrum
            if noutputs==0,
                for ireal=1:nreal,
                    x=double(z(:,1,ireal));
                    switch Order
                        case 1
                            [pxx,f]= pwelch (x,  windFunction, nOverlap, nFFT, fs);
                        case 2
                            error ('high order auto-spectra not defined yet')
                    end
                    if ireal==1,
                        domega=f(2)-f(1);
                        S.dataSet=(pxx);
                        nme = z.chanNames;
                        set (S,'dataSet',pxx,'domainIncr',domega,  'chanNames',{'Gxx'}, ...
                            'comment',['Spectrum of ' nme{1}],'domainName','Frequency', ...
                            'nFFT',nFFT, 'spectType','AUTO');
                    else
                        S.dataSet=cat(3,S.dataSet,pxx);
                        
                    end
                end
                
            elseif ninputs ==1 & noutputs ==1,
                a=double(z);
                x=a(:,1);
                y=a(:,2);
                switch Order
                    case 1
                        [pxy,f]=cpsd(x, y, windFunction, nOverlap, nFFT, fs);
                    case 2
                        error ('high order spectra not supported yet')
                end
                domega=f(2)-f(1);
                set (S,'dataSet',pxy);
                set (S,'domainIncr',domega, 'domainName','Frequency', ...
                    'chanNames',{'Pxy'}, 'nFFT',nFFT,'spectType','Cross');
                
                
            else
                error ('multiple outputs not yet implement');
            end
          
        end
        
        function A = powArea(Sin,varargin)
            % Estimates the power are for determined bands, it can handle
            % multiple bands at the same time as long as the parameters
            % have the same size
            optionList= { { 'fmin' [] 'lower limit of the band'}, ...
                { 'fmax' [] 'upper limit of the band'}...
                {'min_flag' [] 'include the lower limit in area'},...
                {'max_flag' [] 'include the upper limit in area'}};
            arg_parse(optionList, varargin); 
            
            % Check the size of each parameter
            nLims=[size(fmin(:),1),size(fmax(:),1),size(min_flag(:),1),size(max_flag(:),1)];
            % Make sure that the parameters have the same size
            U=unique(nLims(nLims>0));
            if length(U)>1
                error('All the parameters must have the same size')
            end
            
            p=Sin.dataSet;
            f=Sin.domainValues;
            % Empty parameters are set to default values
            if isempty(fmin)
                fmin=zeros(U,1);
            end
            if isempty(fmax)
                fmax=max(f)*ones(U,1);
            end
            if isempty(min_flag)
                min_flag=false(U,1);
            end
            if isempty(max_flag)
                max_flag=false(U,1);
            end
            A=nldat;
            comment=repmat({'power area in the band '},U,1);
            chan=cell(1,size(Sin,2));
            for i=1:size(Sin,2)
                chan{i}=['A' num2str(i)];
            end
            
            % Includes or not lower limit in area
            idx_band=false(size(f,1),U);
            pow=zeros(U,size(p,2),size(p,3));
            
            % For each defined band
            for i=1:U
                if min_flag(i)
                    idx_band(:,i)=(f>=fmin(i));
                    comment{i}=[comment{i} num2str(fmin(i))  '<=f' ];
                else
                    idx_band(:,i)=(f>fmin(i));
                    comment{i}=[comment{i} num2str(fmin(i))  '<f' ];
                end
                % Includes or not upper limit in area
                if max_flag(i)
                    idx_band(:,i)=idx_band(:,i)&(f<=fmax(i));
                    comment{i}=[comment{i} '<=' num2str(fmax(i))];
                else
                    idx_band(:,i)=idx_band(:,i)&(f<fmax(i));
                    comment{i}=[comment{i} '<' num2str(fmax(i))];
                end
                pow(i,:,:)=trapz(f(idx_band(:,i)),p(idx_band(:,i),:,:));
            end
            
            A.dataSet=pow;
            A.comment=comment;
            A.chanNames=chan;
            A.domainName='Area';
        end
        
    end
    
end
