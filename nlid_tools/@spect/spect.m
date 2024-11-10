classdef spect < nldat
    % spect - power spectrum class for NLID toolbox.
    % returns upper and lower cofidnce bounds if paramter
    % confidenceInterval is not NaN. 

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

            s.parameterSet(5)=param('paramName','confidenceLevel', ...
                'paramDefault',NaN, ...
                'paramHelp','Confidence level of estimate', ...
                'paramType','number'  ...
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
                % Ayrio spectra
                for ireal=1:nreal,
                    x=double(z(:,1,ireal));
                    switch Order
                        case 1
                            if isnan(confidenceLevel)
                                [pxx,f]= pwelch (x,  windFunction, nOverlap, nFFT, fs);
                            else
                                [pxx,f,pxxc]= pwelch (x,  windFunction, nOverlap, nFFT, fs, 'ConfidenceLevel', confidenceLevel);
                                pxx=cat(2,pxx,pxxc);
                            end

                        case 2
                            error ('high order auto-spectra not defined yet')
                    end

                    if ireal==1,
                        if ~isnan(confidenceLevel)

                            chanNames={'Gxx' 'lowerConfidence' 'upperConfidence'};
                        else
                            chanNames={'Gxx'};
                        end
                        domega=f(2)-f(1);
                        S.dataSet=(pxx);
                        nme = z.chanNames;
                        set (S,'dataSet',pxx,'domainIncr',domega,  'chanNames',chanNames, ...
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

        function h=plotConfidence (S, plotType, plotColor)
            % Plot spectrum with upper and lower confidence bounds
            % Plot Type can be linear, semilog, loglog
            %Inputs:  S - spectrum to plot
            % plotType - type of plot (['linear' 'semilog' 'loglog'
            % plotColor - 
            % h= cell arrau of pointers to the three lines. 
            if nargin<2
                plotColor='r';
            end
            if nargin==1 | isempty(plotType)
                plotType='linear';
            end
            assign(S.parameterSet); 
            if isnan(confidenceLevel)
                error ('Confidence levels not defined')
            end
            x=domain(S);
            y=S.dataSet;

            switch lower(plotType)
                case 'linear'
                    y=S.dataSet;
                case 'semilog'
                    y=log10(S.dataSet);
                    S.chanUnits={'log power'};
                case 'loglog'
                    x=log10(x); 
                    y=log10(S.dataSet);
                    S.chanUnits={'log power'};
                    S.domainName= {['log ' S.domainName]};
            end

               

            h1= plot(x,S.dataSet(:,1),'color',plotColor);
            h2=line(x,S.dataSet(:,2)); set(h2,'Color',plotColor,'LineStyle', '--');
            h3=line(x,S.dataSet(:,3)); set(h3,'Color',plotColor,'LineStyle', '--');
            h =  {h1 h2 h3};
            legend(S.chanNames);
            xlabel(S.domainName);
            ylabel(S.chanUnits);
            title([S.comment ': spectrum and ' num2str(confidenceLevel) ' confidence level'])
        end

        function [A, paramSet] = powArea(Sin,varargin)
            % Estimates the power are for determined bands, it can handle
            % multiple bands at the same time as long as the parameters
            % have the same size.
            % Input:
            %       - Sin: input Spect object, it can have mutiple channels
            %       and realizations
            % Output:
            %       - A: nldat object containing power area of each
            %       defined band for each channel and realization.
            %       If multiple bands are defined, the areas are stored in
            %       the same order in the A object as different samples.
            %       The channels and realizations of the A object
            %       correspond to those of the input spect object
            %       - paramSet: structure of parameters defined in the
            %       input or as default.
            % Options:
            %       - fmin: lower limit of the band for area estimation. It
            %       can be a single value or a vector for multiple bands
            %       - fmax: upper limit of the band for area estimation. It
            %       can be a single value or a vector for multiple bands
            %       - min_flag: defines if the lower limit is included or
            %       not in the area estimation. If min_flag is false,
            %       fmin<f; if min_flag is true, fmin<=f.
            %       - max_flag: defines if the upper limit is included or
            %       not in the area estimation. If max_flag is false,
            %       f<fmax; if max_flag is true, f<=fmax.
            % ---------------------------------------------------------
            % Note: If multiple bands are analyzed simultaneously, all the
            % parameters must have the same length
            % ---------------------------------------------------------
            optionList= { { 'fmin' [] 'lower limit of the band'}, ...
                { 'fmax' [] 'upper limit of the band'}...
                {'min_flag' [] 'include the lower limit in area'},...
                {'max_flag' [] 'include the upper limit in area'}};
            arg_parse(optionList, varargin);

            % Check the size of each parameter
            nLims=[size(fmin(:),1),size(fmax(:),1),...
                size(min_flag(:),1),size(max_flag(:),1)];
            % Make sure that the parameters have the same size
            U=unique(nLims(nLims>0));
            if length(U)>1
                error('All the parameters must have the same size')
            end

            % Gets the domain and range of the spect object
            p=Sin.dataSet;
            f=Sin.domainValues;
            % Empty parameters are set to default values
            if isempty(fmin)
                % Default value for the minimum frequency is 0
                fmin=zeros(U,1);
            end
            if isempty(fmax)
                % Default value is the maximum frequency of the spect
                % object
                fmax=max(f)*ones(U,1);
            end
            if isempty(min_flag)
                % Default is not to include the lower limit
                min_flag=false(U,1);
            end
            if isempty(max_flag)
                % Default is not to include the upper limit
                max_flag=false(U,1);
            end

            % Generates nldat object to contain the Areas
            A=nldat;
            % Comment on nldat object will outline the analyzed bands
            comment=repmat({'power area in the band '},U,1);
            % Set channel names
            chan=cell(1,size(Sin,2));
            for i=1:size(Sin,2)
                chan{i}=['A' num2str(i)];
            end

            % Indices that are included in each band for area estimation
            idx_band=false(size(f,1),U);
            % Allocates matrix to contain the area estimations
            pow=zeros(U,size(p,2),size(p,3));
            % Defimes the paramSet structure.
            paramSet.fmin=[];
            paramSet.fmax=[];
            paramSet.min_flag=[];
            paramSet.max_flag=[];
            % Allocates the size of the paramSet structure.
            paramSet=repmat(paramSet,U,1);

            % For each defined band
            for i=1:U
                % Define the samples that will be used in the estimation
                % according to the band limits, and the flags
                if min_flag(i)
                    % Includes the minimum frequency
                    idx_band(:,i)=(f>=fmin(i));
                    comment{i}=[comment{i} num2str(fmin(i))  '<=f' ];
                else
                    % Does not include the minimum frequency
                    idx_band(:,i)=(f>fmin(i));
                    comment{i}=[comment{i} num2str(fmin(i))  '<f' ];
                end

                if max_flag(i)
                    % Includes the maximum frequency
                    idx_band(:,i)=idx_band(:,i)&(f<=fmax(i));
                    comment{i}=[comment{i} '<=' num2str(fmax(i))];
                else
                    % Does not include the maximum frequency
                    idx_band(:,i)=idx_band(:,i)&(f<fmax(i));
                    comment{i}=[comment{i} '<' num2str(fmax(i))];
                end
                % Estimates the power within the band for all the channels
                % and realizations simultaneously
                pow(i,:,:)=trapz(f(idx_band(:,i)),p(idx_band(:,i),:,:));

                % Stores the parameters in the paramSet structure
                paramSet(i).fmin=fmin(i);
                paramSet(i).fmax=fmax(i);
                paramSet(i).min_flag=min_flag(i);
                paramSet(i).max_flag=max_flag(i);
            end

            % sets the estimates areas, comment, and channel names to the A
            % object
            A.dataSet=pow;
            A.comment=comment;
            A.chanNames=chan;
            A.domainName='Area';
        end

    end

end
