classdef polynom < nltop
    % polynom - polynomial class for NLID toolbox.
    % polyRange[-1 1] - inpout range
    %       - used with tcheb to scale input to range of tceb polynomials. This must be
    %         set when creating a polynomial apriori
    % polyMean [0] - mean of polynomial
    % polyStd[1] - standard deviation of polynomial
    %         - polyMean and polySTD are used with hermite polynomials to scsle the input to
    %          have zero mean and unti variance. These values must be set to the value expected when
    %          when creating hemrite polynomials a priori.
    % polyOrderSelectMode
    %   'auto' - polynomial order selected on the basis of the minimum data
    %   length
    %   'full' - polynomial order is set to polyOrderMax
    %   'manual' - polynomial order is selected interactively.
    
    % Note that Bsplines and laguerre currently only support linearly
    % incresing domain values
    properties
        polyCoef= nan;
        nInputs=1;
        polyOrder=5;
        polyRange=[-1;1];
        polyMean=0;
        polyStd=1;
        polyVaf=nan;
        parameterSet=param;
    end
    
    methods
        function p = polynom (pin,varargin)
            % add object  specific parameters
            p.parameterSet(1) =param('paramName','polyType','paramDefault','hermite', ...
                'paramHelp','polynomial type',...
                'paramType','select','paramLimits',{'hermite','power','tcheb' 'Bspline' 'laguerre' 'interp1'});
            p.parameterSet(2) =param('paramName','polyOrderSelectMode','paramDefault','auto',...
                'paramHelp','order selection method',...
                'paramType','select','paramLimits',{'auto','full','manual'});
            p.parameterSet(3) =param('paramName','polyOrderMax','paramDefault',10,...
                'paramHelp','maximum order to evaluate');
            polyType=get(p,'polyType');
            
            p.comment='polynomial model';
            if nargin==0;
                return
            elseif isa(pin,'polynom')
                p=nlmkobj(pin,varargin{:});
            else
                p=nlmkobj(p,pin, varargin{:});
            end
            
        end
        
        function sys = set (sys, varargin)
            % Overlaid set function for polynom objects
            v=varargin;
            if iscell(varargin{1}) & length(varargin{1})>1,
                varargin=varargin{1};
            end
            if length(varargin)==1,
                varargin = { varargin{1} []};
            end
            for i=1:2:length(varargin),
                Prop=varargin{i};
                Value=varargin{i+1};
                % Check to see if it is a name is a property and if so set
                % it to Value
                if ismember(Prop, properties(sys)),
                    sys.(Prop)=Value;
                    if strcmp(Prop,'polyCoef'),
                        sys.polyOrder=length(Value)-1;
                    end
                elseif ismember('parameterSet',fieldnames(sys))
                    % Check to see if it is a parameter value
                    % Must change so that handles parameters of different
                    % names.
                    ps = sys.parameterSet;
                    outPs=setval(ps,Prop,Value);
                    if ~isempty(Value) & strcmp(Prop,'polyType'),
                        switch Value
                            case 'Bspline'
                                % Bsplines require addtional parameters to
                                % define them
                                outPs(4) =param('paramName','splineCenters','paramDefault',[0:.2:1]',...
                                    'paramHelp','Centers of splines');
                                outPs(5) =param('paramName','splineSD','paramDefault',.1,...
                                    'paramHelp','Standard Deviation of Splines');
                                sys.polyRange=[0 5];
                                sys.polyCoef= ones(6,1);
                                sys.polyOrder=6;
                                outPs=setval(outPs,'polyOrderMax',sys.polyOrder);
                                outPs=setval(outPs,'polyOrderSelectMode','full');
                            case 'laguerre'
                                % Laguerre needs and addtional paramter to
                                % define it
                                outPs(4) =param('paramName','alfa','paramDefault',.5,...
                                    'paramHelp','alpha parameter for Lagurre polynomials');
                                outPS=outPs(1:4);
                            case 'interp1'
                                outPs=outPs(1);
                                outPs(2) =param('paramName','interpMethod','paramDefault','linear', ...
                                'paramHelp','Interpolation method',...
                                 'paramType','select','paramLimits',{'linear','cubic','spline'});
                                outPs=outPs(1:2);
                            otherwise
                                outPs=outPs(1:3);
                        end
                        
                    elseif ~isempty(Value) & strcmp(Prop,'splineCenters'),
                        %sys.polyCoef=Value*0+1;
                        sys.polyOrder=length(Value);
                        sys.polyRange = [ min(Value) max(Value)];
                    end
                    sys.parameterSet=outPs;
                end
            end
            if ~isempty(inputname(1)),
                assignin('caller',inputname(1),sys);
            end
        end
        
        function p= set.polyCoef (p, value)
            try
                assign(p.parameterSet);
                if strcmp(polyType,'Bspline'),
                    % There is a mismatch between # of coefficient and spline
                    % definitions.
                    if length(value)~=length(splineCenters)
                        error ('Number of coefficients does not match spline definition');
                    end
                end
            catch
            end
            if ndims(value)>2 | ~isreal(value) ,
                error('coefficients must be a real numbers.')
            end
            p.polyCoef = value;
            if p.nInputs == 1;
                p.polyOrder = length(value)-1;
            elseif length(value)==1,
                p.polyOrder=0;
            else
                order = 0;
                ncoeff = 1;
                while order <= 20
                    order = order + 1;
                    ncoeff = ncoeff*(p.nInputs+order)/order;
                    if length(value)==ncoeff
                        p.polyOrder = order;
                        break;
                    end
                end
                if order == 21
                    error('coefficient vector length is incorrect');
                end
            end
        end
        
        
        function V= basisfunction (P, x)
            % Returns basis functions for polynominal P evaluated over x
            if nargin==1,
                pRange=P.polyRange;
                x=(pRange(1):.01:pRange(2))';
            end
            nSamp=length(x);
            assign(P.parameterSet)
            polyOrder=P.polyOrder;
            polyType=get(P,'polyType');
            switch polyType
                case 'hermite'
                    v=multi_herm(x,polyOrder);
                case 'power'
                    v=multi_pwr (x,polyOrder);
                case 'tcheb'
                    v=multi_tcheb(x,polyOrder);
                case 'Bspline'
                    v=generate_B_splines ( x, splineCenters, splineSD );
                case 'laguerre'
                    v=generate_laguerre_basis(nSamp, polyOrder, alfa);
            end
            V=nldat(v,'comment',['Basis functions for ' polyType 'polynomials'],'domainValues',x,'domainName','Input value');
            [nsamp,nchan]=size(V);
            for i=1:nchan,
                chanName{i}=[' Basis ' num2str(i-1)];
            end
            V.chanNames=chanName;
            
            
            
        end
        
        function mdot = ddx(m);
            
            % computes the derivative of the polynomial, m(x), with respect to its
            % argument.
            %
            % syntax:  mdot = ddx(m);
            %
            % where m and mdot are polynom objects
            mdot = m;
            assign(m.parameterSet);
            % turn the polynomial into a power series
            mdot = nlident(mdot,'polyType','power');
            coeff = mdot.polyCoef;
            order = length(coeff);
            % differentiate the power series term by term
            if order==1,
                coeff=0;
            else
            coeff = coeff(2:order).*[1:order-1]';
            end
            set(mdot,'polyCoef',coeff);
            % change back to the original polynomial type.
            mdot = nlident(mdot,'polyType',polyType);
        end
        
        function y = double(x)
            y = cat(1, x.polyRange, x.polyCoef);
        end
        
        function h = hessian(p,z,varargin);
            % computes hessian for least squares estimate of a polynomial
            % make sure that the second argumment is either nldat or double,
            % and contains enough colums for input(s) and output
            
            assign(p.parameterSet);
            if ~(isa (z,'double') | isa(z,'nldat'))
                error (' second argument must be of class double, or nldat');
            end
            [nsamp,nchan,nreal]=size(z);
            nin =p.nInputs;
            
            order = p.polyOrder;
            
            dz=double(z);
            if nchan == nin,
                x=double(domain(z));
                y=dz(:,1:nin);
            else
                x=dz(:,1:nin);
                y=dz(:,nin+1);
            end
            switch lower(polyType);
                case 'power'
                    [w,f]=multi_pwr(x,order);
                case 'hermite'
                    % scale inputs to 0 mean, unit variance.
                    for i=1:nchan-1,
                        x(:,i) = (x(:,i) - polyMean(x(:,i)))/polyStd(x(:,i));
                    end
                    [w,f]=multi_herm(x,order);
                case 'tcheb'
                    % scale inputs to +1 -1;
                    for i=1:nchan-1,
                        a=min(x(:,i));
                        b=max(x(:,i));
                        cmean=(a+b)/2;
                        drange=(b-a)/2;
                        x(:,i)=(x(:,i) - cmean)/drange;
                    end
                    [w,f]=multi_tcheb(x,order);
                otherwise
                    error(['hessian not supported for polyType:' polyType]);
            end
            
            
            h = (w'*w)/nsamp;
            
        end
        function j = jacobian(p,z,varargin);
            % computes jacobian for least squares estimate of a polynomial
            % make sure that the second argumment is either nldat or double,
            % and contains enough colums for input(s) and output
            
            assign(p.parameterSet);
            if ~(isa (z,'double') | isa(z,'nldat'))
                error (' second argument must be of class double, or nldat');
            end
            [nsamp,nchan,nreal]=size(z);
            nin = p.nInputs;
            
            order = p.polyOrder;
            dz=double(z);
            x=dz(:,1:nin);
            switch lower(polyType);
                case 'power'
                    [j,f]=multi_pwr(x,order);
                case 'hermite'
                    % scale inputs to 0 mean, unit variance.
                    for i=1:nchan-1,
                        x(:,i) = (x(:,i) - polyMean(x(:,i)))/polyStd(x(:,i));
                    end
                    [j,f]=multi_herm(x,order);
                case 'tcheb'
                    % scale inputs to +1 -1;
                    for i=1:nchan-1,
                        a=min(x(:,i));
                        b=max(x(:,i));
                        cmean=(a+b)/2;
                        drange=(b-a)/2;
                        x(:,i)=(x(:,i) - cmean)/drange;
                    end
                    [j,f]=multi_tcheb(x,order);
                otherwise
                    error(['jacobian not support for polyType:' polyType])
                    
            end
        end
        
        function y = nlsim ( sys, xIn )
            % polynom/nlsim - simulate response of polynominal series to input data set
            inputClass=class(xIn);
            nin=sys.nInputs;
            assign(sys.parameterSet);
            [~,nchan,nreal]=size(xIn);
            segdat_flag = 0;
            switch inputClass
                case'segdat'
                    segdat_flag = 1;
                    y=xIn;
                case 'nldat'
                    y=xIn;
                case 'double'
                    xIn=nldat(xIn);
                    y=xIn;
            end
            if (nchan ~= nin)
                error ('number of input channels does not match polynominal');
            end
            x=xIn.dataSet;
            [nr,nc]=size(x);
            if nin > nc,
                error ('dimension mismatch');
            else
                xData=xIn.dataSet;
                for iReal=1:nreal,
                    x=xData(:,:,iReal);
                    
                    switch lower(polyType)
                        
                        case 'power'
                            p=multi_pwr(x,sys.polyOrder);
                            yout=p*sys.polyCoef(:);
                        case 'hermite'
                            for i=1:nchan,
                                x(:,i) = (x(:,i) - sys.polyMean(i))/sys.polyStd(i);
                            end
                            p=multi_herm (x,sys.polyOrder);
                            yout=p*sys.polyCoef(:);
                        case 'tcheb'
                            for i=1:nchan,
                                % scale input with same scale factor
                                % used in estimating the polynomial
                                r=sys.polyRange;
                                a=r(1,i);
                                b=r(2,i);
                                x(:,i)=(2*x(:,i) -(b+a))/(b-a);
                            end
                            p=multi_tcheb (x,sys.polyOrder);
                            yout=p*sys.polyCoef(:);
                        case 'bspline'
                            p=generate_B_splines(x,splineCenters, splineSD);
                            yout=p*sys.polyCoef;
                        case 'laguerre'
                            p=generate_laguerre_basis(length(x),sys.polyOrder, alfa);
                            yout=p*sys.polyCoef;
                        case 'interp1'
                            yout=interp1(sys.polyCoef(:,1),sys.polyCoef(:,2),x,interpMethod);
                    end
                    Y(:,:,iReal)=yout;
                    
                end
                
            end
           
            set(y,'dataSet',Y,'comment',[ polyType ' series prediction'], 'domainIncr', ...
                xIn.domainIncr,'domainStart',xIn.domainStart);
            % polynom/nlsim
        end
        %%
        function nlmtst (p)
            polynomDemo
        end
        function plot (p, z)
            % overload plot function for polynom objects
            % z - optional input showing orginal poly set
            nin=p.nInputs;
            range=p.polyRange;
            if isnan(range);
                umean = p.polyMean;
                ustd = p.polyStd;
                if isnan(umean+ustd)
                    range=[-1 1];
                else
                    range=[umean-3*ustd umean+4*ustd];
                end
            end
            
            if nin==1,
                x=linspace(range(1),range(2))';
                y=double(nlsim(p,x));
                plot (x,y);
                if nargin > 1,
                    hold on
                    
                    plot (z(:,1).dataSet, z(:,2).dataSet,'r.');
                    hold off
                end
                
                title(p.comment);
                xlabel('input');
                ylabel('output');
            elseif nin==2,
                x=linspace(range(1,1), range(2,1))';
                y=linspace(range(1,2), range(2,2))';
                [x,y]=meshgrid(x,y);
                z=cat(2,x(:),y(:));
                zp=nlsim(sys,z);
                z=double(zp);
                z=reshape (z,length(x),length(y));
                mesh (x,y,z);
            else
                error('plotting not available for polynominal with > 2 inputs');
            end
        end
        
        function zOut = rdivide(z,x)
            % Divide coefficients of a polynomial by scalzr
            zOut=z;
            zOut.polyCoef=z.polyCoef./x;
        end
        function zOut = times(z,x),
            % Multiply the coefficients of a polynomial by a scalar
            zOut=z;
            zOut.polyCoef=z.polyCoef.*x;
        end
        
        
        
        
        
        function p = nlident (pin, z, varargin );
            % polynom/nlident - Overlaid nlident for polynom class
            
            global POLY_ORDER POLY_DONE
            
            if isa(z,'char')
                % z and varargin are property value pairs
                p = nltransform(pin,{z ,varargin{1:end}});
                return
            end
            
            
            if ~(isa (z,'double') | isa(z,'nldat'))
                error (' Second argument must be of class double, or nldat');
            end
            p=pin;
            if nargin >2
                set(p,varargin{:});
            end
            assign (p.parameterSet);
            
            [nsamp,nchan,nreal]=size(z);
            nin = p.nInputs;
            dz=double(z);
            if nchan == nin,
                x=double(domain(z));
                y=dz(:,1:nin);
            else
                x=dz(:,1:nin);
                y=dz(:,nin+1);
            end
            
            
            set (p,'polyRange', [min(x) ; max(x)]);
            set (p,'polyMean',mean(x),'polyStd',std(x));
            switch lower(polyType);
                case 'power'
                    [W,f]=multi_pwr(x,polyOrderMax);
                case 'hermite'
                    % Scale inputs to 0 mean, unit variance.
                    for i=1:nchan-1,
                        x(:,i) = (x(:,i) - mean(x(:,i)))/std(x(:,i));
                    end
                    [W,f]=multi_herm(x,polyOrderMax);
                case 'tcheb'
                    % Scale inputs to +1 -1;
                    for i=1:nchan-1,
                        a=min(x(:,i));
                        b=max(x(:,i));
                        cmean=(a+b)/2;
                        drange=(b-a)/2;
                        x(:,i)=(x(:,i) - cmean)/drange;
                    end
                    [W,f]=multi_tcheb(x,polyOrderMax);
                case 'bspline'
                    % set defaults center and SD if  not defined
                    t=double(z(:,1));
                    xMin=min(t);
                    xMax=max(t);
                    p.polyRange(1)=xMin;
                    p.polyRange(2)=xMax;
                    deltx=(xMax-xMin)/(p.polyOrder-1)
                    splineCenters=[xMin:deltx:xMax]';
                    splineSD=deltx;
                    set(p,'splineCenters',splineCenters);
                    set(p,'splineSD',splineSD);
                    polyOrderMax=length(splineCenters);
                    set(p,'polyOrderMax',polyOrderMax);
                    bf=generate_B_splines(t,splineCenters,splineSD);
                    W=double(bf);
                case 'laguerre'
                    t=double(z(:,1));
                    xMin=min(t);
                    xMax=max(t);
                    p.polyRange(1)=xMin;
                    p.polyRange(2)=xMax;
                    disp('laguerre polynomials fit to ramp time data');
                    L=generate_laguerre_basis(length(x), polyOrderMax,alfa);
                    W=double(L);
                    polyOrderSelectMode='full'
                case 'interp1'
                    polyOrder=pin.polyOrder;
                    x=double(z(:,1));
                    [xs,is]=sort(x);
                    y=double(z(:,2));
                    ys=y(is);
                    xMin=unique(min(xs));
                    xMax=unique(max(xs));
                    deltaX=(xMax-xMin)/polyOrder;
                    curVal=xMin;
                    coef=zeros(polyOrder,2);
                    for i=1:polyOrder+1
                        j=min(find(ge(xs,curVal)))
                        curVal=curVal+deltaX;
                        coef(i,1)=xs(j);
                        coef(i,2)=ys(j);
                    end
                    set(p,'polyCoef',coef,'polyRange',[xMin xMax]);
                    return
            end
            [Q,R] = qr(W,0);
            QtY = Q'*y;
            y_sum_sq = sum(y.^2);
            
            
            % Choose polynomial order....
            
            
            switch polyOrderSelectMode
                case 'auto'
                    [coefs,mdls] = poly_mdls(y_sum_sq,nsamp,R,QtY,polyOrderMax,nin);
                    [val,pos] = min(mdls);
                    order = pos-1;
                    d = factorial(nin+order)/(factorial(nin)*factorial(order));
                    coef = coefs(1:d,order+1);
                case 'manual'
                    [coefs,mdls,vafs] = poly_mdls(y_sum_sq,nsamp,R,QtY,polyOrderMax,nin);
                    poly_gui('init',mdls,vafs);
                    done = 0;
                    while ~done;
                        done = POLY_DONE;
                        pause(0.25);
                    end
                    order = POLY_ORDER;
                    d = factorial(nin+order)/(factorial(nin)*factorial(order));
                    coef = coefs(1:d,order+1);
                case 'full'
                    % use maximum order,
                    order = polyOrderMax;
                    npar = length(QtY);
                    coef = qr_solve(R,QtY,npar);
                otherwise
                    error('unrecognized mode');
            end
            
            %% B spine
            
            if strcmp(polyType,'Bspline'),
                set (p,'polyCoef',coef);
                set (p,'polyOrder',order);
            else
                set (p,'polyCoef',coef,'polyOrder',order);
                
            end
            
            return
            
        end
        
        function [outSNL, info] = snl2hwr(SNL,u,varargin) 
            %% Apr. 12th, 2023: This function converts a polynomial Static Nonlinearity (SNL) model into a Half-Wave Rectifier (HWR) SNL model.
            %% The input SNL is a polynom object with any of these polyTypes: 'hermite','power','tcheb' 'Bspline' 'laguerre' 'interp1'
            %% The output SNL is a polynom object with polyType equal to 'interp1'. We represent a HWR in polynom object as an 'interp1' look up table.
            % This function uses a grid search method to fit an optimal HWR to a polynomial SNL, 
            % where optimality can be defined using various configurable options:
            %    + optimization/minimization criterion defined with argument 'error_critetion' that can take values: 'curve_err' and 'pred_err'
            %    + minimization error metrics defined with argument 'error_metric that can take values: 'rmse' (root mean saquared error) and 'mae' (max absolute error)
            
            %-- Inputs: 
            %           SNL: An polynom object from NLID library 
            %             u: An nldat object containing the input to the SNL used for its identification 
            %                  
            %-- Outputs:
            %           outSNL: Best HWR fit HWR captured into a polynom object with polyType equal to 'interp1' 
            %           info: A structure containing information/metadata about the conversion including:
            %                 .err: A status text message indicating whether HWR fit worked or not
            %                 .optimErr: A field to plot error surface on the grid of slope and threshold; it has 3 subfields:
            %                            .values, .meshThreshod, .meshSlope, 
            %                            
            
            options={{'slope_range' [] 'the range of slopes to search for'} ...
                     {'n_slopes' 150 'number of slopes to try in the slope range'} ...
                     {'n_thresholds' 150 'number of thresholds to try in the threshold range'} ...
                     {'optim_criterion' 'pred_err' 'Minimization error criterion. Available options are: pred_err (SNL output prediction error), curve_err (SNL curve estimation error)'} ...
                     {'error_metric' 'rmse' 'Minimization error metrics. Available options are: rmse (root mean squared error), mae (max absolute error)'} ...
                     {'input_prctiles' [1, 99] 'the input range, in percentile, to use for threshold search and curve_err minimization'} ...
                     {'plot_figures' false 'the flag to plot or suppress plotting the figures. Default is zero, which means suppress or not plot.'} ...
                 };
            
            if arg_parse(options,varargin)
                return
            end
            
            if isempty(slope_range)
                info.err = 'You must specify a range for slopes to search for.';
                info.optimErr = [];
                outSNL = [];
                disp(info.err)
                return
            end
            
            %% Determine whether the SNL is a polymial object i 
            if ~isa(SNL,'polynom')
                info.err = 'The input SNL is not a polynomial object of NLID.';
                info.optimErr = [];
                outSNL = [];
                disp(info.err)
                return
            end
            
            %% Find the threshold and slope range to search 
            %  threshold range is always found from input data range defined by input_prctiles  
            min_t = prctile(u.dataSet,input_prctiles(1));
            max_t = prctile(u.dataSet,input_prctiles(2));
            
            %% Main body of the function
            %% Creating grids for Threshold and Slope 
            t = linspace(min_t,max_t,n_thresholds);
            
            max_a = max(slope_range); 
            min_a = min(slope_range);
            a = linspace(min_a,max_a,n_slopes);
            
            %% Loop over slope and threshold to calculate the error metric based on the selected error criterion
            %% First, we need to calculate the response of the SNL to a ramp input and the input used for identification
            nRampPoints = 500; 
            uRamp = nldat(linspace(min_t,max_t,nRampPoints)');
            yRamp = nlsim(SNL,uRamp);        %-- Response of the SNL to a ramp input
            
            uData = u.dataSet;
            
            switch optim_criterion
                case 'pred_err'
                    %-- With this optimization criterion, the optimization output must be NL response to u
                    y = nlsim(SNL,u);  %-- Response of the SNL to u
                    y = y.dataSet;
                
                case 'curve_err'
                    %-- With this optimization criterion, the optimization output must be NL response to ramp 
                    y = yRamp.dataSet; 
            
                otherwise
                    info.err = 'The provided optimization criterion is not supported.';
                    info.optimErr = [];
                    outSNL = [];
                    disp(info.err);
                    return;
            end
            
            %-- To simplify notation, set x to be equal to the ramp in the input range
            x = uRamp.dataSet;
            xMax = max(x);
            xMin = min(x);
            
            errMetric = zeros(n_slopes,n_thresholds);
            
            for i = 1:length(a) 
                for j = 1:length(t)
                    yhat = zeros(size(y));
                    smallerThanT = x<=t(j);
                    switch optim_criterion
                        case 'curve_err'                   
                            %-- Find the baseline
                            baseline = mean(y(smallerThanT));
                            
                            %-- Calculate the HWR curve
                            yhat(smallerThanT) = baseline;
                            biggerThanT = x>t(j);
                            yhat(biggerThanT) = a(i)*(x(biggerThanT)-t(j)) + baseline;
                        
                        case 'pred_err' 
                            %-- Find the baseline
                            yr = yRamp.dataSet;
                            baseline = mean(yr(smallerThanT));
                            
                            %-- Simulate the response of the HWR (captured as a polynom object) to input u
                            %-- The reason I did not convert this 
                            for k = 1:length(uData)
                                if uData(k) <= t(j)
                                    yhat(k,1) = baseline;
                                else
                                    yhat(k,1) = a(i)*(uData(k)-t(j)) + baseline;
                                end
                            end
            
                        otherwise
                            info.err = 'This error criterion is not implemented. Available error criteria are: pred_err and curve_err';
                            info.optimErr = [];
                            outSNL = [];
                            disp(info.err);
                            return;
                    end
                    
                    switch error_metric
                        case 'rmse'
                            errMetric(i,j) = sqrt(mean((y-yhat).*(y-yhat)));
                        case 'mae'
                            errMetric(i,j) = max(abs((y-yhat).*(y-yhat)));
                        otherwise
                            info.err = 'This error metric is not implemented. Available error metric are: rmse and mae';
                            info.optimErr = [];
                            outSNL = [];
                            disp(info.err)
                            return;
                    end
            
                end
            end
            
            %% Finding the best half-wave rectifier fit (slope, threshold, baseline)
            [minRowVal,minRowInd] = min(errMetric);
            [~,minColInd] = min(minRowVal);
            
            slope = a(minRowInd(minColInd)); 
            threshold = t(minColInd); 
            
            smallerThanT = find(x<=threshold);
            switch optim_criterion
                case 'curve_err'
                    baseline = mean(y(smallerThanT));
                case 'pred_err'
                    baseline = mean(yr(smallerThanT));
                otherwise
                    info.err = 'The provided optimization criterion is not defined.';
                    info.optimErr = [];
                    outSNL = [];
                    disp(info.err);
                    return;
            end
            
            %% Setting the output SNL
            outSNL = polynom('polyType','interp1');
            epsilon = 1e-12;   %-- The reason for adding/subtracting this very small epsilon is to avoid an error with 'interp1' as it cannot receive two identical coeffients.
            if threshold == xMin 
                threshold = xMin + epsilon;
            elseif threshold == xMax
                threshold = xMax - epsilon;
            end
            
            polyCoeff_X = [xMin, threshold, xMax]';
            polyCoeff_Y = [baseline, baseline, slope*(xMax-threshold) + baseline]';
            
            set(outSNL,'polyCoef',[polyCoeff_X,polyCoeff_Y]);
            
            %== Setting the output info structure
            [meshT,meshS] = meshgrid(t',a');
            info.err = 'No errors! The HWR was successfully estimated.';
            info.optimErr.values = errMetric;
            info.optimErr.meshThreshold = meshT;
            info.optimErr.meshSlope = meshS;
            
            if plot_figures
                %== Plot the figure comparing the polynomial SNL with its HWR fit
                yHWR = nlsim(outSNL,uRamp);
                ySNL = nlsim(SNL,uRamp);
                figure; plot(x,ySNL.dataSet,'b',x,yHWR.dataSet,'r'); xlabel('UA'); ylabel('NL Output'); legend('Ident SNL','HWR Fit'); title(['Optimization criterion was: ',optim_criterion])
            end
        
        end    


    end
end
%%

function B = generate_B_splines(q,centers,sd)
% q = domain
% centers = centers for splines
% sd - standard deviations for each spline
q=q(:);
centers=centers(:);
sd=sd(:);
number=length(centers);
B=zeros(length(q),number);
for i=1:number
    B(:,i)=(1/(2*pi*sd)^(1/2))*exp(-(0.5/(sd^2))*((q-centers(i)).*(q-centers(i))));
end

return

end


%%
function b = generate_laguerre_basis (irf_len,max_order,alfa)
%% This function generates the Laguerre orthonormal basis functions
%++ Author: Ehsan Sobhani (10 April 2014)
%++ This is based on Marmaleris book OR formula (11) of his paper titled:
%++ "Identification of Nonlinear Biological Systems Using Laguerre Expansions of Kernels", Annals of Biomed. Eng., vol. 21, pp. 573-589, 1993.

%++ The inputs are:
% 1) irf_len (integer number of samples)
% 2) max_order (integer). The maximum order of Laguerre expansion.
% 3) alfa (0<alfa<1). Laguerre parameter.

%++ The outputs are:
% 1) b (real with size [irf_len,max_order+1])

%=== Initialization
b = zeros(irf_len,max_order+1);

for t = 1:irf_len
    for j = 0:max_order
        gain = alfa^((t-j)/2) * (1-alfa)^(0.5);
        summation = 0;
        for k = 0:j
            argument = (-1)^k * combination(t,k) * combination(j,k) * alfa^(j-k) * (1-alfa)^k;
            summation = summation + argument;
        end
        b(t,j+1) = gain * summation;
    end
end

end
%%
function p = nltransform(pin,inputs)
% Transform polynomial from one basis to another

ni = length(inputs);

if rem(ni,2)~=0,
    error('Property/value pairs must come in even number.')
end


p = pin;

assign (pin.parameterSet);
for i = 1:2:ni,
    % Set each PV pair in turn
    Property=inputs{i};
    Value = inputs{i+1};
    set(p, Property,Value);
    if  ~isnan(p.polyCoef)
        switch Property
            case 'polyOrder'
                % change order and pad/truncate coeff
                
                NumInputs = p.nInputs;
                OldOrder = pin.polyOrder;
                OldCoeffs = pin.polyCoef;
                NewCoeffs = multi_pwr(zeros(NumInputs,1),Value);
                if NewOrder > OldOrder
                    OldLength = length(OldCoeffs);
                    NewCoeffs(1:OldLength) = OldCoeffs;
                else
                    NewLength = length(NewCoeffs);
                    NewCoeffs = OldCoeffs(1:NewLength);
                end
                polyCoef = NewCoeffs;
                
            case 'PolyType'
                % change type and recompute coeff
                
                if strcmp(PolyType,Value)
                    % Type is already set so do nothing
                    break
                elseif p.nInputs > 1,
                    error ('Conversion not yet implemented for multiple inputs')
                else
                    Pstats = [p.polyRange(2) p.polyRange(1) p.polyMean p.polyStd];
                    Pstats = AllStats(Pstats);
                    cold=(piN.Coef);
                    OldType = PolyType;
                    cnew = poly_convert(cold, Pstats, OldType,Value);
                    p.polyCoef=cnew;
                end
                
            otherwise
                % set behaves appropriately, so use it.
                
        end
    end
end
end


function PStats = AllStats(PStats)
% fills in missing input stats

if isnan(PStats(1))
    % range has not yet been specified
    % use 3 sigma bounds around mean (if set)
    if isnan(PStats(3))
        % mean not set either, assume zero mean
        PStats(3) = 0;
    end
    if isnan(PStats(4))
        PStats(4) = 1;
    end
    PStats(1) = PStats(3)+3*PStats(4);
    PStats(2) = PStats(3)-3*PStats(4);
end

if isnan(PStats(3))
    % mean has not been specified
    % let the max and min be 3 sigma around the mean.
    PStats(3) = (PStats(1)+PStats(2))/2;
    PStats(4) = (PStats(1)-PStats(3))/3;
end

return
end


function [coefs,mdls,vafs] = poly_mdls(yss,LenY,R,QtY,OrderMax,nin);

mdls = zeros(OrderMax+1,1);
vafs = mdls;


max_coefs =  factorial(nin+OrderMax)/(factorial(nin)*factorial(OrderMax));
coefs = zeros(max_coefs,OrderMax+1);


% MDL penalty gain
k = log(LenY)/LenY;

if nin == 1
    ends = 1+[0:OrderMax]';
else
    ends = zeros(OrderMax+1,1);
    % ends go at (n+q)! / n! q!
    ends(1) = 1;
    for q = 1:OrderMax
        ends(q+1) = ends(q) * (nin+q)/q;
    end
end

RTR = R'*R;

% compute coefs, vafs and mdls
for q = 0:OrderMax
    d = ends(q+1);
    coef = qr_solve(R,QtY,d);
    coefs(1:d,q+1) = coef;
    outvar = coef'*RTR(1:d,1:d)*coef;
    resid = yss-outvar;
    mdls(q+1) = resid*(1+k*d);
    vafs(q+1) = 100*(outvar/yss);
end


return
end


function coef = qr_solve(R,QtY,npar)
% solves for npar coefficients using precomputed QR factorization
% intened to be used to estimate reduced order models
% R is the R matrix from a QR decomposition
% QtY is Q^T y, where Q is from the QR decomposition, and y is the output.

QtY = QtY(1:npar);
R = R(1:npar,1:npar);
coef = R\QtY;

end

% copyright 2003, robert e kearney
% this file is part of the nlid toolbox, and is released under the gnu
% general public license for details, see ../copying.txt and ../gpl.txt