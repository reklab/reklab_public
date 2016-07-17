function system= hammer_subspace (z,ps)
% This function estimates a Hammerstein between input and output stored in columns of z.
% This function uses MOESP subspace approach
% 
% syntax:  system= hammer_subspace (z)
% where:
%       hammer is a nlbl object containing a hammerstein cascade.
%       z is an nldat objects containing the input-output data
%
%
% The available options (and their default values) are:
%          'nLagLE' (16):                  length of linear component IRF in samples
%          'maxOrderNLE' (12):             maximum order for nonlinearity
%          'nHankle' (20):                 Size of hankle matrix; must be larger than the linear system order
%          'threshSSE:  Threshold on SSE for terminitiaon of the iterative search
%          'orderSelectMethod' ('manual'): Method for order selection; can be manual or auto
%          'nDelayInput' (0):              Input delay (samples)
%
% See
% [1] Jalaleddini, K., and Kearney, R.E., Subspace Identification of SISO Hammerstein Systems:
% Application to Stretch Reflex Identification, 2013, in IEEE Transactions on Biomedical Engineering (under review).
%                       &&
% [2] Jalaleddini, K., and Kearney R.E., An Iterative Algorithm for the Subspace Identification
% of SISO Hammerstein Systems, in Proceedings of 18th World Congress IFAC,
% pp.11779-11784, 2011.

% options={ ...
%          {'nLagLE' 16 'length of linear component IRF in samples'} ...
%          {'maxOrderNLE' 12 'maximum order for nonlinearity'} ...
%          {'nHankle' 20 'Size of hankle matrix; must be larger than the linear system order'} ...
%          {'threshSSE ' 10^-10 'Threshold on % change of SSE for terminitiaon of the iterative search'} ...
%          {'orderSelectMethod' 'manual' 'Method for order selection; can be manual or auto'} ...
%          {'nDelayInput' 0 'Input delay (samples)'} ...
%      };
% if arg_parse_c('exact',options,varargin');
%      return
% end
%% Initialization
assign(ps);
it = 1000; % Maximum number of iterations 
%z = z - mean(z);
ts = get(z,'domainIncr');
input = z(:,1);
input = input.dataSet;
input = del(input,nDelayInput);
output = z(:,2);
%Normalizing input for Tchebychev expansion
avg = (max(input)+min(input))/2;
rng = max(input) - min(input);
un = (input - avg)*2/rng;
%u_r is the Tchebychev expansion of input
u_r = multi_tcheb(un,maxOrderNLE-1);
%% Estimate the extended observability matrix for u_r as the input and output
[S, R] = dordpi(u_r,output.dataSet,hankleSize);
%Selecting the order of the linear system
if strcmp('preset',orderSelectMethodLE),
    n = orderLE;
else
n = orderselect(S,orderSelectMethodLE);
end
%% Estimate the A and C matrix from extended observability matrix
[AT, CT] = destac(R,n);
%Return if the system poles are outside of the unit circle
if ~isempty(find(abs(eig(AT))>1, 1))
    disp('Identified System is unstable. Function returns void matrices...')
    static_nl = polynom;
    system = nlbl;
    system_ss = ssm;
    set(system,'elements',{static_nl,system_ss},'idMethod','subspace');
    return
elseif n==0
    disp('Zero-order system. Function returns void matrices...')
    static_nl = polynom;
    system = nlbl;
    system_ss = ssm;
    set(system,'elements',{static_nl,system_ss},'idMethod','subspace');
    return
else
%% Estimate Phi using dcalcphi from SMI toolbox
%Phi is the regressor of equation (14-15) of [1]. Phi is the regressor for
%the overparameterized model
[Phi,T] = dcalcphi(u_r,output.dataSet,AT,CT);
%Remove the last column of row since it has a copy of output
Phi = Phi(:,1:size(Phi,2)-1);
%% Estimate BT and DT elements and nonlinearity coefficients using the iterative algorithm
%alpha contains the coefficient of the Tchebychev nonlinearity
%bd contains the elements of the B and D matrices
alpha = ones(maxOrderNLE,1);
alpha = alpha/norm(alpha);
s1 = 1;
s2 = 1;
DS=output.dataSet;
sse=DS'*DS;
%Phi_alpha is the right hand side of (20) in [1].
%Phi_bd is the right hand side of (23) in [1].
%Phi_alpha_it = Phi*Phi_alpha defines the regressor for estimation of B and D elements
%Ph_alpha_it is the left hand side of (20)
%Phi_bd_it = Phi*Phi_bd defines the regressor for estimation of alpha elements
%Ph_bd_it is the left hand side of (23)
%Now performing the iteration. Each iteration has 3 parts: 1)update alpha;
%2) update bd; 3) normalization
for i=1:it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the regressor Phi_alpha_it based on alpha for estimation update of bd
    Phi_alpha = zeros((n+1)*maxOrderNLE,n+1);
    for j = 1 : maxOrderNLE
        Phi_alpha ((j-1)*n+1:j*n,1:n) = eye(n) * alpha(j);
    end
    Phi_alpha (maxOrderNLE*n+1:end,n+1) = alpha;
    Phi_alpha_it = Phi * Phi_alpha;
%Use least-sqaures to update bd and compute the corresponding SSE
    bd = Phi_alpha_it \ DS;
    sse_c = DS'*DS - bd'*Phi_alpha_it'*DS;
    nsse_c=sse_c/sse;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the regressor Phi_bd_it based on bd to update the estimate of of alpha
    Phi_bd=zeros((n+1)*maxOrderNLE,maxOrderNLE);
    for j=1:maxOrderNLE
        Phi_bd((j-1)*n+1:j*n,j) = bd(1:n);
    end
    Phi_bd (maxOrderNLE*n+1:end,:) = eye(maxOrderNLE)*bd(n+1);
    Phi_bd_it = Phi*Phi_bd;
%Use least-sqaures to update alpha and compute the corresponding SSE    
    alpha = Phi_bd_it \ DS;
    sse_b = DS'*DS - alpha'*Phi_bd_it'*DS;
    nsse_b = sse_b/sse;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize bd and alpha for the next round of iteration
    h = sign(alpha(1));
    bd = h*bd * norm(alpha);
    alpha = alpha / norm(alpha)*h;
%Terminate the iteration if the improvement in NSSE is less than threshSSE
    if displayFlag,
        disp([s1 s1-nsse_c s2 s2-nsse_b]);
    end
    if (s1-nsse_c<threshNSE) && (s2-nsse_b<threshNSE)
        break
    end
    
    s1 = nsse_c;
    s2 = nsse_b;
    
end
%The B matrix has the first n elements of bd
BT = bd(1:n);
BT = T*BT;
%The D matrix has the last element of bd
DT = bd(n+1);
DT = DT';

%% Construct the output nlbl object
%Construct the estimated nonlinearity based on the Tchebychev coefficient
newMin = min(input);
newMax = max(input);
newMean = mean(input);
newStd = std(input);
alpha = polynom('polyCoef',alpha,'polyType','tcheb',...
    'comment','Static Nonlinearity','polyRange',[newMin;newMax],...
    'polyMean',newMean,'polyStd',newStd);
%Constrcut ssm object for the linear system
system_ss = ssm;
set(system_ss,'A',AT,'B',BT,'C',CT,'D',DT,'domainIncr',ts,'nDelayInput',nDelayInput);
%Concatenate the polynom and IRF objects to construct the estimated Hammerstein system
static_nl = alpha;
system = nlbl;
set(system,'comment','NL System identified using subspace method ', ...
    'elements',{static_nl,system_ss},'idMethod','subspace', ...
    'orderLE', n);
end
end
function [d_x] = del(x,nDelay)

% this function adds a delay to the input signal 
% x is the input signal
% nDelay - dealy in samples
% d_x is the delayed singal
% n is the length of delay

n = length(x);
d_x = zeros(size(x));
for i=1:n-nDelay
    d_x(nDelay+i) = x(i);
end
  end
