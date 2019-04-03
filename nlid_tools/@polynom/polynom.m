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
                'paramType','select','paramLimits',{'hermite','power','tcheb' 'Bspline' 'laguerre'});
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
                            otherwise
                                outPs=outPs(1:3);
                        end
                        
                    elseif ~isempty(Value) & strcmp(Prop,'splineCenters'),
                        sys.polyCoef=Value*0+1;
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
            V=nldat(v,'comment',['Basis functions for ' polyType ],'domainValues',x,'domainName','Input value');
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
            coeff = coeff(2:order).*[1:order-1]';
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
                    
            end
        end
        
        function y = nlsim ( sys, xin )
            % polynom/nlsim - simulate response of polynominal series to input data set
            
            
            nin=sys.nInputs;
            assign(sys.parameterSet);
            [~,nchan,nreal]=size(xin);
            segdat_flag = 0;
            if strcmp(class(xin),'segdat')
                segdat_flag = 1;
                xIN = xin;
                xin =nldat(xin(:,1));
                if nchan == 2
                    nchan = 1;
                end
                if (nchan ~= nin)
                    error ('number of input channels does not match polynominal');
                end
                %nchan = nchan - 1;
            elseif (nchan ~= nin)
                error ('number of input channels does not match polynominal');
            end
            if strcmp(class(xin),'nldat')
                y=xin;
            elseif isa(class(xin),'segdat')
                y=segdat;
            else
                y = nldat;
            end
            x=double (xin);
            [nr,nc]=size(x);
            if nin > nc,
                error ('dimension mismatch');
            else
                for iReal=1:nreal,
                    x=double(xin(:,:,iReal));
                    
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
                    end
                    Y(:,:,iReal)=yout;
                    
                end
                
            end
            if segdat_flag
                dataset = xIN.dataSet;
                dataset (:,1) = Y;
                y = xIN;
                set(y,'dataSet',dataset);
            else
                set(y,'dataSet',Y);
            end
            set(y,'comment',[ polyType ' series prediction'],'domainValues',double(xin));
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
                    t=domain(z);
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
                    t=domain(z);
                    xMin=min(t);
                    xMax=max(t); 
                    p.polyRange(1)=xMin;
                    p.polyRange(2)=xMax;
                    disp('laguerre polynomials fit to ramp time data');
                    L=generate_laguerre_basis(length(x), polyOrderMax,alfa);
                    W=double(L);
                    polyOrderSelectMode='full'
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