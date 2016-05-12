function Mout = nlident (Min, z, varargin );
% Overlaid nlident for tvm models 

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


if nargin < 2,
    disp('NLIDtakes two inputs for irf objects: irf, Z' );
elseif nargin > 2,
    set(Min,varargin)
end
assign (Min.Parameters);
Mout=Min;
ModelType= get(Min, 'Model_Type');
%
% Identify time-varying polynomial
%
switch ModelType
case 'polynom'
    [nsamp, nchan, nreal]=size(z);
    disp('TVPolynom');
    [nsamp,nchan,nreal]=size(z);
    p=polynom;
    for i= 1:nsamp;
        x1=squeeze(z(i,1,:));
        y1=squeeze(z(i,2,:)); 
        z1 = cat(2,x1,y1);
        p= polynom(p,z1);
        P{i}=p;
    end
    set (Mout,'Data',P);
    
    %
    % TV IRF
    %
di=get(z,'DomainIncr');set(Mout,'domainincr',di);
dn=get(z,'DomainName');set(Mout,'domainname',dn);
ds=get(z,'DomainStart');set(Mout,'domainstart',ds);
dv=get(z,'DomainValues');set(Mout,'domainvalues',dv);

case 'irf'
    %
    %  Parmaeter "Method" determines the method used:
    %   'tvfil' - finds least-squares solution using the data itself.
    %   'corr' - finds least-squares solution using correlation functions.
    %   'pseudo' - correlation function approach as with 'corr', but
    %   uses Dave's technique (adapted to tv case) to select singular vectors
    %   Default is 'corr'.
    [nsamp,nchan,nr]=size(z);
    P=get(Min,'Parameters');
    assign(P);
    if isnan (NLags)
        NLags=min(64, round(nsamp/100));
    end
    Ts=get(z,'DomainIncr');
    if isnan(Ts),
        error('IRF estimation requires domainincr to be specified');
    end
    x = squeeze(double(z(:,1,:)));
    y =squeeze(double(z(:,2,:)));
    conf_level=NaN;
    [Hident,bound,sing_vectors,cpu] = tv_ident(x,y,Ts,NSides,NLags,Method,conf_level);
    if NSides ==2,
        LagStart=-(NLags-1)*Ts/2;
        TimeStart=0;
    else
        LagStart=0;
        TimeStart=0;
        %
        % end
        
    end
    [nsamp,nlag]=size(Hident);
    I=irf;
    set(I, 'Nsides','NSides','domainincr',Ts,'domainstart',LagStart);
    IZero=I;
    set(IZero,'data',Hident(1,:)'*0);
    j=1;
    for i=1:(NLags-1)/NSides,
        IA{j,1}=IZero;
        j=j+1;
    end
     for i=1:nsamp,
        set(I,'data',Hident(i,:)');
        IA{j,1}=I;
        j=j+1;
    end
    if NSides ==2,
        for i=1:(NLags-1)/NSides,
        IA{j,1}=IZero;
        j=j+1;
    end
    end

    set(Mout,'Model_type','irf','domainincr',Ts,'domainname','time','nsides', NSides, 'nlags', NLags, ...
        'data',IA, 'domainstart',TimeStart,'Method',Method);
 
case 'nlbl'
    % Identify time varying Hammerstein Model
    [nsamp,nchan,nr]=size(z);
    P=get(Min,'Parameters');
    assign(P);
    if isnan (NLags)
        NLags=min(64, round(nsamp/100));
    end
    Ts=get(z,'DomainIncr');
    if isnan(Ts),
        error('IRF estimation requires domainincr to be specified');
    end
    x = squeeze(double(z(:,1,:)));
    y =squeeze(double(z(:,2,:)));
    if NSides ==1,
        M1=0;
        M2=NLags;
    else
        M1=floor(NLags/2);
        M2=M1;
    end
    
  [coeffs,range,IRF,num_iterations] = tv_hammer_i(x,y,Ts,OrderMax,-M1,M2,Tolerance);
if NSides ==2,
        LagStart=-(NLags)*Ts/2;
        TimeStart=NLags*Ts;
        NTop=NLags/2;
        NBottom=NLags/2;
    else
        LagStart=0;
        TimeStart=(NLags)*Ts;
        NTop=NLags;
        NBottom=0;
        %
        % end
        
    end
    p=polynom;
    set(p,'polyType','power','OrderMax',OrderMax,'Order',OrderMax);
    iout=irf;
    set(iout,'NSides',NSides, 'domainincr',Ts, ...
        'NLags',NLags,'domainstart',LagStart, ...
        'domainname','Lag','Comment','IRF');
    NL=nlbl;
    IRF=IRF';
    for i=1:length(coeffs),
        c=coeffs{i};
        c=(fliplr(c))';
      set (p,'Coef',c, 'Range', range(i,:));
      set(iout,'data',IRF(:,i));
      NL{1,1}=p;
      NL{1,2}=iout;
      el{i,1}=NL;
  end
  set(Mout,'domainincr',Ts,'domainname','time', ...
        'data',el, 'domainstart',TimeStart);
 
  Mout=pad(Mout,NTop,NBottom);
  
    %
otherwise
    error ('tvm - nlident bad model type' ) 
    
    
end



