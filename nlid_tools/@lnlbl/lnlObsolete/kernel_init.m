function m0 = kernel_init(m0,uy);
% initial lnl model, based on kernel estimates.
% system must have significant even polynomial terms.

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

%% Moved to lmlbl 
u = uy(:,1);

%1 extract objects
blocks = get(m0,'elements');
h = blocks{1};
m = blocks{2};
g = blocks{3};

disp('kernel init start')
get(m,'polyOrderMax')




Ts = get(u,'domainIncr');

N = length(double(u));
hlen = length(h);
glen = length(g);
memory = hlen+glen-1;

% estimate volterra series model
vs = vseries;
set(vs,'idMethod','foa','nLags',memory);
vs = nlident(vs,uy);
kernels = get(vs,'elements');
K2 = double(kernels{3});


% find first non-zero row in second kernel

% project second-order kernel onto significant singular vectors from input
% auto-correlation function.
phi = cor;
set(phi,'nLags',memory,'corType','covar','biasMode','biased');
phi = nlident(phi,u);
Phiuu = toeplitz(double(phi));
[U,S,V] = svd(Phiuu);
ss = diag(S);
gaps = ss(1:memory-1)./ss(2:memory);
[val,num_svs] = max(gaps)
K2p =  U(:,1:num_svs)*V(:,1:num_svs)'*K2;

K2sumsq = sum(K2p.^2);
thresh = 0.1*max(K2sumsq);
slice = min(find(K2sumsq>thresh))

set(h,'dataSet',K2p(1:hlen,slice)*Ts,'domainIncr',Ts);
set(g,'domainIncr',Ts,'nLags',glen);
set(m0,'elements',{h m g});


disp('kernel init before call to BestHammer')
get(m,'polyOrderMax')



m0 = BestHammer(m0,uy,'hk');
% change this to 'init', or 'ng' once implemented
