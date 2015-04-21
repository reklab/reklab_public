
function [hammer,  xpoly,  xlinear] = tv_hammer_system(duration,dt,verbose)

% hammer = hammer_system(duration,dt,verbose)
%
% Creates tv_hammerstein object used to test identification technique
%
% duration:      desired duration of the time-varying behaviour (in seconds)
% dt:            sampling interval 
%
% hammer:        tv_hammerstein object 
%

% Check that number of input arguments is valid.
num_args = nargin;
error(nargchk(2,3,num_args));

% Check that input arguments are valid.
if ~(duration >0)
   error('duration must be a positive real number')
end
if ~(dt>0)
   error('dt must be a positive real number')
end
if (num_args > 2) % user specified verbose
   if ~(isnumeric(verbose))
      error('verbose must be a numeric entry')
   end
else % user did not specify verbose
   % set default value for verbose
   verbose = 0;
end  % if

% Set constants
h_duration = 0.5;               % duration of dynamics = 500 ms

% Set other parameter values
hlength = ceil(h_duration/dt)+1;  % number of points in each IRF
num_instants = ceil(duration/dt);

% Generate parameters for static nonlinearity
pn=polynom;
set(pn,'type','power');
c1_1 = 0.5; c2_1 = 0.75; c3_1 = 1; c4_1 = 0;
c1_2 = 1.5; c2_2 = 1.25; c3_2 = 1; c4_2 = 0;
s1_1 = 0;   s2_1 = 2;    s3_1 = 2; s4_1 = 0;
s1_2 = 0;   s2_2 = 4;    s3_2 = 2; s4_2 = 0;
v = (0:(num_instants/4-1))';
c1 = c1_1+(c1_2-c1_1)/(num_instants/4-1)*v;
c2 = c2_1+(c2_2-c2_1)/(num_instants/4-1)*v;
c3 = c3_1+(c3_2-c3_1)/(num_instants/4-1)*v;
c4 = c4_1+(c4_2-c4_1)/(num_instants/4-1)*v;
s1 = s1_1+(s1_2-s1_1)/(num_instants/4-1)*v;
s2 = s2_1+(s2_2-s2_1)/(num_instants/4-1)*v;
s3 = s3_1+(s3_2-s3_1)/(num_instants/4-1)*v;
s4 = s4_1+(s4_2-s4_1)/(num_instants/4-1)*v;
% p1 = [c1;s1;flipud(s1);flipud(c1)];
% p2 = [c2;s2;flipud(s2);flipud(c2)];
% p3 = [c3;s3;flipud(s3);flipud(c3)];
% p4 = [c4;s4;flipud(s4);flipud(c4)];
 p1 = [c1;s1;flipud(s1);flipud(c1)];
 p2 = [c2;s2;flipud(s2);flipud(c2)];
 p3 = [c3;s3;flipud(s3);flipud(c3)];
 p4 = [c4;s4;flipud(s4);flipud(c4)];

p = [p1 p2 p3 p4];
for i = 1:size(p,1)
    set (pn,'coef',flipud(p(i,:)'),'range', [-3 3],'Order',length(p(i,:))-1);
    ptv{i,1}=pn;
end  % for i
% Create tv_static_nonlin object 
SN = tvm;
set(SN,'Model_type','polynom','domainincr',dt,'Data',ptv);
xpoly=SN;
%SN = tv_static_nonlin(coeffs,range,'dt',dt,'type','poly');

% Generate IRF's 
i1 = 0.0015; b1 = 0.0216; k1 = 0.3719; % corresponds to 100th row of IRF produced by reflex_pathway.m
i2 = 7.621e-4; b2 = 0.0144; k2 = 0.1249; % corresponds to 200th row of IRF produced by reflex_pathway.m
I = [i1*ones(num_instants/2,1);i2*ones(num_instants/4,1);i1*ones(num_instants/4,1)];
B = [b1*ones(num_instants/2,1);b2*ones(num_instants/4,1);b1*ones(num_instants/4,1)];
K = [k1*ones(num_instants/2,1);k2*ones(num_instants/4,1);k1*ones(num_instants/4,1)];
DL = tv_2(I,B,K,hlength,dt,0.001);
% normalize s.t. area under curve is 1
DLD=get(DL,'data');
for i=1:length(DLD),
    Di=DLD{i};
    d=get(Di,'data');
    d=d/(sum(d)*dt);
    set(Di,'data',d);
    DLD{i,1}=Di;
end
set(DL,'data',DLD,'domainincr',dt);
xlinear=DL;


% Create tv_hammerstein object

NL=nlbl;
for i=1:length(DL),
        xpd=xpoly.data;
        xld=xlinear.data;
    
        set(NL,'Elements', { xpd{i} xld{i} });
 EL{i,1}=NL;
 

end
hammer=tvm;
set(hammer,'model_type','nlbl','data',EL, 'domainincr',dt);
set(hammer,'Nlags', hlength, 'Method','Simulation');
% If verbose, plot some information.
if (verbose)

   figure
   subplot(2,2,1)
   plot(p1)
   title('p1')
   subplot(2,2,2)
   plot(p2)
   title('p2')
   subplot(2,2,3)
   plot(p3)
   title('p3')
   subplot(2,2,4)
   plot(p4)
   title('p4')
   
   figure
   subplot(3,1,1)
   plot(I)
   title('I')
   subplot(3,1,2)
   plot(B)
   title('B')
   subplot(3,1,3)
   plot(K)
   title('K')
   
   figure
   plot(SN)
   title('Static nonlinearity')

   figure
   plot(IRF)
   title('Time-varying dynamics')

end  % if verbose




