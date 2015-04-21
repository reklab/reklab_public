function pc_model = lnl2pc(lnl,hlen,numpaths,method);
% converts lnl into parallel cascde.
%
% syntax: pc_model = lnl2pc(lnl,hlen,numpaths,method);
%
% method can be 'kernels', 'kernel2','delays'

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% 1. construct IRF matrix
% 2. compute polynomials
% 3. put elements into a pcascade object.

if nargin < 4
  method = 'kernels';
end

H = construct_irfs(lnl,numpaths,hlen,method);
M = construct_polys(lnl,H);

pc_model = make_pcascade(H,M);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pc_model = make_pcascade(H,M);
% constructs an object of type pcascade from the cell array of IRF objects,
% H, and the cell array of POLYNOM objects, M. 

% get properties from IRFs
Ts = get(H{1},'domainincr');
numpaths = length(H);
hlen = get(H{1},'Nlags');

% get properties from POLYNOMs
order = get(M{1},'Order');

% initialize the parallel cascade object.
pc_model = pcascade;
set(pc_model,'OrderMax',order,'NPaths',numpaths,'NLags',hlen);
for i = 1:numpaths
  pc_elements{i,1} = H{i};
  pc_elements{i,2} = M{i};
end
set(pc_model,'Elements',pc_elements);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = construct_irfs(lnl,numpaths,hlen,method);
% constructs a cell array containing IRF objects which will be used to form
% a parallel cascade model of the LNL cascde object lnl.

H = {};
Irfs = zeros(hlen,numpaths);
if nargin < 3
  method = 'kernels'
end

% get properties from LNL object
elements = get(lnl,'Elements');
g = get(elements{1},'data');
h = get(elements{3},'data');
Ts = get(elements{1},'domainincr');

% create a matrix containing the irfs, using the given method.
switch method
  case 'kernels'
    k1 = get(lnl2volt(lnl,hlen,1),'data');
    k2 = get(lnl2volt(lnl,hlen,2),'data');
    k1 = k1/sqrt(sum(k1.^2));
    [U,S,V] = svd(k2);
    Irfs = [k1 U(:,1:numpaths-1)]/Ts;
  case 'delays'
    hk = zeros(hlen,1);
    glen = length(g);
    if glen <= hlen
      hk(1:glen) = g;
    else
      hk = g(1:hlen);
    end
    hk = hk/(Ts*sqrt(sum(hk.^2)));
    for i = 1:numpaths
      Irfs(:,i) = hk;
      hk = [0;hk(1:hlen-1)];
    end
  case 'kernel2'
    k2 = get(lnl2volt(lnl,hlen,2),'data');
    [U,S,V] = svd(k2);
    Irfs = U(:,1:numpaths)/Ts;    
end

% create a dummy IRF object, and give it the properties common to all of the
% IRFs. 
hk = irf;
set(hk,'domainincr',Ts,'NLags',hlen);

for i = 1:numpaths
  set(hk,'data',Irfs(:,i));
  H{i} = hk;
end
H = H';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = construct_polys(lnl,H);
% constructs a cell array containing the polynomials associated with the
% IRFs contained in the cell-array H, that best fit the LNL cascade lnl.



% get properties from LNL object
elements = get(lnl,'Elements');
g = get(elements{1},'data');
Ts = get(elements{1},'domainincr');
glen = length(g);
m = elements{2};
mc = nlident(m,'type','power');
Coeffp = get(mc,'Coef');
order = get(mc,'Order');
h = get(elements{3},'data');
hlen = length(h);

numpaths = length(H);
NLags = get(H{1},'NLags');
Irfs = zeros(NLags,numpaths);
for i = 1:numpaths
  Irfs(:,i) = get(H{i},'data');
end



% the n'th order polynomial coefficients will be fitted using
%
%  c = K^{-n}G^nh
%  where
%  K(i,j) = \sum_t Irfs(t,i)Irfs(t,j) 
%  G(i,j) = \sum_t g(i-j)Irfs(t,i)
K = Irfs'*Irfs;
Gdel = zeros(hlen,NLags);
gk = zeros(NLags,1);
if NLags >= glen
  gk(1:glen) = g;
else
  gk = g(1:NLags);
end
for i = 1:hlen
  Gdel(:,i) = gk;
  gk = [0; gk(1:NLags-1)];
end
G = Irfs'*Gdel;

AllCoeffs = zeros(order+1,numpaths);
AllCoeffs(1,1) = Coeffp(1)*Ts*sum(h);
for i = 1:order
  Kn = K.^i;
  Gn = G.^i;
  c = inv(Kn)*Gn*h*Ts;
  AllCoeffs(i+1,:) = Coeffp(i+1)*c';
end

% Now, figure out the input attributes for each of the paths.
AllRanges = zeros(2,numpaths);
AllStds = zeros(numpaths,1);
AllMeans = AllStds;
% 1. assume a white gaussian input the the original LNL cascade, compute the
%    mean and variance that would reproduce the stats for its nonlinearity.
xRange = get(m,'Range');
xStd = get(m,'Std');
% Alternatively, assume Range corresponds to 3 sigma bounds, and compute an
% alternate Std, based on that.
xAltStd = (xRange(2)-xRange(1))/6;
% Be conservative, choose the smaller of the two.
xStd = min(xStd,xAltStd);
xGainW = Ts * sqrt(sum(g.^2));
uStd = xStd/xGainW;


xMean = get(m,'Mean');
xGainDC = Ts * sum(g);
% uMean = xMean/xGainDC;  % we have to make an exception for high pass
% systems. 
uMean = 0;


% 2. Compute the mean and variance of the outputs of each of the linear
%    elements, given this input.  Use these as the input properties for the
%    corresponding nonlinearity.

for i = 1:numpaths
  hk = Irfs(:,i);
  DCGain = Ts * sum(hk);
  Gain = Ts * sqrt(sum(hk.^2));
  wMean = uMean*DCGain;
  wStd = uStd*Gain;
  wMin = wMean - 3*wStd;
  wMax = wMean + 3*wStd;
  AllRanges(:,i) = [wMin; wMax];
  AllStds(i) = wStd;
  AllMeans(i) = wMean;
end

% And finally, assemble all of the polynomials into the output cell-array. 

mk = polynom;
set(mk,'Type','power','Order',order);
for i = 1:numpaths
  set(mk,'Coef',AllCoeffs(:,i),'Range',AllRanges(:,i),...
      'Std',AllStds(i),'Mean',AllMeans(i));
  M{i} = mk;
end
M = M';






