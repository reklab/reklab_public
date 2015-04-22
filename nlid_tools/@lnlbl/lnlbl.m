classdef lnlbl < nlm
    % lnlbl - linear-nonlinear-linear block model class for NLID toolbox.
     
    properties
    end
    
    methods
        function LNL = lnlbl  (a,varargin)
            
            LNL.parameterSet(1)=param('paramName','idMethod', ...
                'paramDefault','hk', ...
                'paramHelp',...
                'Iterative method used to refine initial estimate',...
                'paramType' ,'select','paramLimits', {'hk', 'init'});
            LNL.parameterSet(2)= param('paramName','initMethod',...
                'paramDefault','order1', ...
                'paramHelp', 'Initial identification method', ...
                'paramType','select',...
                'paramLimits',{'kernels','self','order1'});
            
            LNL.parameterSet(3)=param('paramName','nLags1', ...
                'paramDefault',NaN,...
                'paramHelp',  'Number of lags in first linear element' ,...
                'paramType','number', ...
                'paramLimits', {0 1000});
            
            LNL.parameterSet(4)=param('paramName','polyOrderMax', ...
                'paramDefault',2,...
                'paramHelp','Maximum order for polynomisl nonlinearity' ,...
                'paramType','number', ...
                'paramLimits', {0 10});
            
            LNL.parameterSet(5)=param('paramName','nLags2', ...
                'paramDefault',NaN,...
                'paramHelp','Number of lags in second linear element' ,...
                'paramType','number', ...
                'paramLimits', {0 1000});
            
            LNL.parameterSet(6)=param('paramName','pseudoInvMode', ...
                'paramDefault','full', ...
                'paramHelp', 'pseudo-inverse order selection mode ', ...
                'paramType','select',...
                'paramLimits',{'full','auto','manual'});
            
            LNL.parameterSet(7) = param('paramName','hkTolerance', ...
                'paramDefault',0.01, ...
                'paramHelp','Improvement required to continue iteration', ...
                'paramType','number',...
                'paramLimits', {0 100});
            
            LNL.parameterSet(8) = param('paramName','nhkMaxIts', ...
                'paramDefault',20, ...
                'paramHelp','Maximum number of iterations', ...
                'paramType',   'number',...
                'paramLimits',{1 1000});
            
            LNL.parameterSet(9) = param('paramName','nhkMaxInner', ...
                'paramDefault',10, ...
                'paramHelp','Maximum number of iterations between Hammerstein updates',...
                'paramType','number', ...
                'paramLimits',{1 100});
            
            LNL.parameterSet(10)=param('paramName','hkAccel', ...
                'paramDefault',0.8, ...
                'paramHelp','ridge multiplied by decel after unsuccessful update',...
                'paramType','number', ...
                'paramLimits', {0.001 0.999});
            
            LNL.parameterSet(11)=param('paramName','hkDecel', ...
                'paramDefault',2, ...
                'paramHelp','ridge size multipled by accel after successful update',...
                'paramType','number', ...
                'paramLimits', {1.0001 inf});
            
            LNL.parameterSet(12)=param('paramName','hkInitialStep', ...
                'paramDefault',10,'paramHelp',...
                'initial stepsize','paramType','number',...
                'paramLimits',{0 inf});
            i=irf;
            j=irf;
            t=polynom;
            set(t,'polyType','tcheb');
            elements = { i t j};
            set (LNL,'elements',elements);
            if nargin==0;
                return
            elseif nargin==1,
                LNL=nlmkobj(LNL,a);
            elseif isa(a,'lnlbl')
                LNL=nlmkobj(a,varargin{:});
            else
                NLN=nlmkobj(NLN ,a,varargin{:});
            end
        end
        
        function bl  = nlident (bl, z, varargin)
            % Identify a lnlbl
            if nargin > 2
                set (bl,varargin);
            end
            if isa(z,'nldat') | isa(z,'double')
                % if numlags is undefined, for either block, set a sensible default.
                numlags=get(bl,'nLags1');
                if isnan(numlags'),
                    numlags= max(32,length(z)/100);
                    set(bl,'nLags1',numlags);
                end
                numlags=get(bl,'nLags2');
                if isnan(numlags'),
                    numlags= max(32,length(z)/100);
                    set(bl,'nLags2',numlags);
                end
                if isa(z,'nldat')
                    Ts=get(z,'domainIncr');
                else
                    subsys = get(bl,'elements');
                    f1 = subsys{1};
                    Ts = get(f1,'domainIncr');
                    z = nldat(z,'domainIncr',Ts);
                end
                subsys = get(bl,'elements');
                g = subsys{1};  % First IRF
                p = subsys{2};  % Polynomial
                h = subsys{3};  % Second IRF
                
                glen = get(bl,'nLags1');
                hlen = get(bl,'nLags2');
                mode = get(bl,'pseudoInvMode');
                set(g,'nLags',glen,'irfPseudoInvMode',mode);
                set(h,'nLags',hlen,'irfPseudoInvMode',mode);
                order = get(bl,'polyOrderMax');
                set(p,'polyOrderMax', order);
                set(bl,'elements',{g p h});
                
                
                x=z(:,1);
                y=z(:,2);
                switch lower(get(bl,'initMethod'))
                    case 'kernels'
                        bl = kernel_init(bl,z);
                    case 'self'
                        % do nothing
                    case 'order1'
                        bl = hk_init(bl,z);
                    otherwise
                        error('Unsupported Initialization Method');
                end
                
                
                switch lower(get(bl,'idMethod'))
                    case 'init'
                        set(bl,'comment','Initial LNL Model estimate');
                    case 'hk'
                        bl = hk_ident(bl,z);
                    otherwise
                        error('Unsupported Refinement Method');
                end
                
            else
                error('conversions to models of class lnbl not yet implemented');
            end
            
        end
        
        function hs = smo (h, npass);
            % smo for nlnbl objects
            hs=h;
            el=get(h,'elements');
            irf1=smo(el{1},npass);
            irf2=smo(el{3},npass);
            el{1}=irf1;
            el{2}=irf2;
            set (hs,'elements',el);
        end
        
        
    end
end

function m0 = hk_init(m0,uy);
% initial lnl model, first-order fit to first-order kernel
% for details see
%
%      M.J. Korenberg and I.W. Hunter, The identification of nonlinear
%      biological systems: LNL cascade models.
%      Biological Cybernetics, 55:125-134, 1986.
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt


%1 extract objects
blocks = get(m0,'elements');
h = blocks{1};
m = blocks{2};
g = blocks{3};

hlen = get(m0,'nLags1');
glen = get(m0,'nLags2');
memory = hlen+glen;
set(h,'nLags',hlen);
set(g,'nLags',glen);

% estimate first order kernel
k1 = h;
set(k1,'nLags',memory,'irfPseudoInvMode','auto');
k1 = nlident(k1,uy);
Ts = get(k1,'domainIncr');

% fit an appropriate first-order unity gain system to the kernel.
% use 10/90 rise time in step response as about 2.2 time constants.
% using the 10/90 rise time should eliminate any problems caused by
% delays in the system. (However, it won't work for high-pass systems)

step_response = cumsum(double(k1));
fin_value = step_response(end);
if fin_value<0
    ste_response = -step_response;
    fin_value = -fin_value;
end

t90 = min(find(step_response > 0.9*fin_value));
t10 = min(find(step_response > 0.1*fin_value));
tau = (Ts/2.2)*(t90-t10);

t = [0:hlen-1]'*Ts;
if tau > 0
    hh = exp(-t/tau)/Ts;
else
    hh = [1/Ts;zeros(hlen-1,1)];
end
set(h,'dataSet',hh,'domainIncr',Ts);
set(g,'domainIncr',Ts);
set(m0,'elements',{h m g});
m0 = BestHammer(m0,uy,'hk');
end

function [m0,err,vf] = BestHammer (m0,uy,method);
% computes Hammerstein model to follow first IRF.

if nargin < 3
    idMethod = 'hk';
end
blocks = get(m0,'elements');
h = blocks{1};
m = blocks{2};
g = blocks{3};

u = uy(:,1);
y = uy(:,2);

Ts = get(u,'domainIncr');
order = get(m,'polyOrderMax');
hlen = get(g,'nLags');
x = nlsim(h,u);

hdata = cat(2,x,y);
mh = nlbl(hdata,'idMethod',method,'polyOrderMax',order,'nLags',hlen,...
    'iterationTolerance',0.01);
hblocks = get(mh,'elements');
m = hblocks{1};
g = hblocks{2};
set(m0,'elements',{h m g});

if nargout > 1
    testout = nlsim(m0,u);
    err = y - testout;
    vf = double(vaf(y,testout));
end
end

function m0 = hk_ident(m0,uy);
% Identification of LNL cascade using the Hunter-Korenberg
% iterative cross-correlation method.
%
% See:  M.J. Korenberg and I.W. Hunter, The identification of nonlinear
%      biological systems: LNL cascade models.
%      Biological Cybernetics, 55:125-134, 1986.
%
% syntax  m0 = hk_ident(m0,uy);
%   TP = [num_svs, accel, decel, NumIts1,NumIts2];
%
% num_svs: number of svs used in computation of initial IRF
%               if num_svs == 0, then start from m0


% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt

% get data
u = uy(:,1);
y = uy(:,2);


%1 extract objects
blocks = get(m0,'elements');
h = blocks{1};
m = blocks{2};
g = blocks{3};
accel = get(m0,'hkAccel');
decel = get(m0,'hkDecel');
NumIts2 = get(m0,'nhkMaxInner');
NumIts1 = get(m0,'nhkMaxIts');
step_size = get(m0,'hkInitialStep');
mh = nlbl;
set(mh,'elements',{m g});
Ts = get(u,'domainIncr');
N = length(double(u));
hlen = length(double(h));

testout = nlsim(m0,u);
err = y - testout;
vf = double(vaf(y,testout));
old_vaf = vf;
vafs = vf;
new_h = h;
new_m = m0;


for i = 1:NumIts1
    a_count = 0;
    for j = 1:NumIts2
        stuff = get(m0,'elements');
        h = stuff{1}; m = stuff{2}; g = stuff{3};
        delta = LMStep(u,err,m0,step_size);
        set(new_h,'dataSet', double(h) + delta);
        set(new_m,'elements',{new_h m g});
        testout = nlsim(new_m,u);
        vf = double(vaf(y,testout));
        if vf > old_vaf
            old_vaf = vf;
            err = y - testout;
            m0 = NormalizeLNL(new_m);
            step_size = accel*step_size;
            a_count = a_count + 1;
            %      disp('accelerating');
        else
            step_size = decel*step_size;
            %      disp('decelerating');
        end
    end
    
    disp([num2str(a_count) ' accelerations']);
    if a_count > 0
        [new_m,new_err,vf] = BestHammer (m0,uy,'hk');
        a_count = 0;
        disp(vf-old_vaf)
        if vf > old_vaf
            old_vaf = vf;
            m0 = new_m;
            err = new_err;
            disp('Hammerstein Updated');
            %      step_size = 100*step_size;
        else
            disp('No Improvement in Hammerstein Fit');
        end
    end
    vafs(i+1) = old_vaf;
    disp(old_vaf);
end
end


function [m0,h] = NormalizeLNL(m0);


blocks = get(m0,'elements');
h = blocks{1};
m = blocks{2};
g = blocks{3};

hh = double(h);
gain = std(hh);
set(h,'dataSet',hh/gain);

range = get(m,'polyRange');
set(m,'polyRange',range/gain);

gg = double(g);
gain = std(gg);
set(g,'dataSet',gg/gain);

coeffs = get(m,'polyCoef');
coeffs = coeffs*gain;
set(m,'polyCoef',coeffs);

set(m0,'elements',{h m g});
end






function delta = LMStep(u,err,m0,step_size);

N = length(double(u));
Ts = get(u,'domainIncr');
stuff = get(m0,'elements');
h = stuff{1};
hlen = length(h);

m = stuff{2};
mdot = ddx(m);
xx = double(nlsim(mdot,nlsim(h,u)));
g = stuff{3};
X = zeros(N,hlen);
ud = double(u);
for i = 1:hlen
    temp = nldat(ud.*xx,'domainIncr',Ts);
    X(:,i) = Ts*double(nlsim(g,temp));
    ud = [0;ud(1:N-1)];
end

%delta = inv(X'*X + step_size*eye(hlen))*X'*double(err);
delta = [X;sqrt(step_size)*eye(hlen)]\[double(err);zeros(hlen,1)];
end

function X = FDJacobian(u,m0)

N = length(double(u));
Ts = get(u,'domainIncr');
stuff = get(m0,'elements');
h = stuff{1};
hlen = length(h);

m = stuff{2};
g = stuff{3};

X = zeros(N,hlen);
y = nlsim(m0,u);
yd = y;
md = m0;
hd = h;


for i = 1:hlen
    hh = double(h);
    delta = 0.001*hh(i);
    hh(i) = hh(i) + delta;
    set(hd,'dataSet',hh);
    set(md,'elements',{hd m g});
    yd = nlsim(md,u);
    X(:,i) = double(yd - y)/delta;
end
end








