function [Z,M,  zFull] = nlid_sim (option,U,varargin)
% NLID_SIM - simulate various nonlinear systems
% usage:
%			[z,m,n]  = nlid_sim (option,u, pltflg)
% outputs: 
%  z - simulated data set
%  m - model simulated 
%  zFull - input, noise free output, noise
% inputs: 
%	option 	= system type
%       LP1 - firat order  low pass
%       HP1 - first order high pass
%       L1 - one-sided low pass filter
%       H1 - seocond-order high pass
%       H2 - Bessel high pass
%       L2 - second-order low pass
%       LNRECT - Wiener system - second order low-opass followed by full wave rectifier
%       LN2 - Wiender system with quadratic nonlinearity
%       LN3 - Wiener systems with cubiv nonlinerity
%       N2L - Hammerstein system with quadratic nonlinearity
%       N3L - Hammerstien system with cubic nonlinearity.
%       N3HP - Hammerstein system:  cubic nonlinearity + HP filter
%       PC - parallel cascade model
%       POLY - static nonlineariy
%       STATIC_LINEAR  - static linear system
%       Cuber - cuber 
%   U - input signal 
%   varargin
%       pltflg [true/false]	= plot model simulated
%       'nsmo' '0' ' tnumber of times to smooth [0-N]'} ...
%       'delay_time' 0 'Output delay (sec))'} ...
%       'domain_incr' .01 'Sampling increment (sec))'} ...
%       'noise_level' 0 'Noise std/output STD'
% $Revisions: %
% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt

paramList={{'plotFlg' 'false' 'plot [true/false]'} ...
    {'nsmo' '0' ' tnumber of times to smooth [0-N]'} ...
    {'delay_time' 0 'Output delay (sec))'} ...
    {'domain_incr' .01 'Sampling increment (sec))'} ...
    {'noise_level' 0 'Noise std/Signal STD'}
    };

if arg_parse(paramList,varargin);
    return
end

if nargin<2  || isempty(U)
    u=4*(rand(2200,1)-.5);
    domain_incr=0.01;
    U=nldat(u,'domainIncr',domain_incr);
    U=smo(U,1); 
else
    domain_incr=U.domainIncr;
end
if nargin==0 | isempty(option)
    opts = { 'LP1' 'HP1'  'L1' 'L2'  'LNRect' 'LN2' 'LN3' 'N3L' 'N2L' ...
        'LNL' 'PC' 'POLY' 'Static_Linear' 'Cuber' 'N3HP'};
    o=menu('Option',opts);
    option=opts{o};
end
comment='undefined';
k=1;
z=.5;
w=5*2*pi;
t=0:domain_incr:.5;
a1=[ k z w];
a2=[ 1 .75 25*2*pi];
[irf1]=fkzw(t,a1);
I1=irf;
set(I1,'dataSet',irf1,'domainIncr',domain_incr);
[irf2, di]=fkzw(t,a2);
I2=irf;set(I2,'dataSet',irf2,'domainIncr',domain_incr);
P2=polynom;
set(P2,'polyType','power','polyCoef',[10 25 10 ]');
P3=polynom('polyType','power','polyCoef',[200 50 -10 5]');


%  And the coefficients of the Nonlinearity are...

option=upper(option);
switch option
    % linear
    case {'LP1' 'HP1'}
        t=0:domain_incr:.1;
        freq=15*2*pi;
        TC=1/freq;
        irfLP =domain_incr*freq*exp(-t*freq);
        I=irf;
        switch option
            case 'LP1'
                set(I,'dataSet',irfLP', 'domainIncr',domain_incr);
                Y=nlsim(I,U);
          
            case 'HP1'
                irfHP=-irfLP;
                irfHP(1)=-sum(irfHP);
                set(I,'domainIncr',domain_incr,'dataSet',irfHP');
                Y=nlsim(I,U);
           
        end
        M=I;
        
    case 'L1'
        Y = nlsim(I1,U);
        set (Y,'comment', 'One-side low pass');
        comment='L1';
        M=I1;
    case 'H1'
        [irf1]=fkzw(t,[1 .75 15*2*pi]);
        I1=irf;
        set(I1,'dataSet',irf1,'domainIncr',domain_incr);
        Y = U - nlsim(I1,U);
        set (Y,'comment', 'One-side low pass');
        comment='L1';
        M=I1;
    case 'H2'
        [B,A]=butter(1,.2,'high');
        Y=filter(U,B,A);
        set(Y,'comment','Bessel filter');
           
    case 'L2'
        Y = nlsim(I2,U);
        M=I2;
        subplot (1,1,1);
        plot(I2);
        title('L2');
        comment='Linear data set 2';
        % LN
        
    case 'LNRECT'
        disp('ln')
        X = nlsim(I1,U);                % filter with the first L
        Y=abs(X);
        comment='LN Threshold data set';
        
    case 'LN2'
        x = nlsim(I1,U);
        Y= nlsim(P2,x);
        M=lnbl; set(M,'elements',{I1 P2});
        comment='LN2 data set';
    case 'LN3'
        x = nlsim(I1,U);             % filter with the first L
        Y= nlsim(P3,x);
        M=lnbl; set(M,'elements',{I1 P3});
        comment='LN3 data set';
        % NL
    case 'N3L'
        x=nlsim(P3,U);
        Y=nlsim(I1,x);
        M=nlbl; set(M,'elements',{P3 I1});
        comment='NL Cubic data set';
    case 'N3HP'
        % Hamerstein system with high pass dynamics
        % Higpass IRf
        u=nldat(randn(1000,1),'domainIncr',domain_incr);
        y=nlsim(I1,u);
        z=u-y;
        irfHigh= irf( cat(2,u,z),'nSides',2,'nLags',51);
        M=nlbl('elements', { P3 irfHigh});
        comment='NL system with highpass dynamics'; 
        Y=nlsim(M,U); 

    case 'N2L'
        x=nlsim(P2,U);
        u=double(U); 
        set(P2,'polyRange',[ min(u) max(u)]);
        Y=nlsim(I1,x);
        M=nlbl; set(M,'elements',{P2 I1});
        comment='NL Quadratic data set';
        %LNL
    case 'LNL'
        x=nlsim(I1,U);
        z=nlsim(P3,x);
        Y=nlsim(I2,z);
        comment='LNL data set';
        M=lnlbl; set(M,'elements',{I1 P3 I2});
    case 'PC'
        x1=nlsim(I1,U);
        y1=nlsim(P2,x1);
        x2=nlsim(I2,U);
        y2=nlsim(P3,x2);
        Y=y1-y2;
        M=pcascade;
        set(M,'elements',{ I1 P2; I2 P3});
        Y=x1+x2;
        comment='Parallel Cascade Data set';
    case 'POLY'
        Y=nlsim(P3,U);
        M=P3;
    case 'STATIC_LINEAR'
        Y=5*U;
        M='Linear';
    case 'CUBER'
        Y=U*U*U
        M='cube';
    otherwise
        error (['Option not defined:' option]);
end
%% Delay data if necessary
if delay_time>0,
    delayN=round(delay_time/domain_incr);
    if delayN > length(U),
        error('Specified delay greater than data length');
    end
    y =Y.dataSet;
    yNew=circshift(Y.dataSet,delayN);
    yNew(1:delayN)=zeros(delayN,1);
    set(Y,'dataSet',yNew);
end
cleanY=Y; 
%% Add noise
Noise=Y*0;
if noise_level>0,
    y=Y.dataSet;
    noise =randn(length(y),1);
    noise=(noise_level*sqrt(var(y)/var(noise)))*noise;
    yNew=y+noise;
    set(Y,'dataSet',yNew);
    set (Noise, 'dataSet',noise);
end
   

Z=cat(2, U,Y);
zFull= cat(2, U, cleanY, Noise);
set(zFull,'chanNames',{'zin' 'yout' 'noise'});




set (Z, 'chanNames', { 'zin' 'zout' } );
set (Z,'comment',comment,'domainName','Time (s)');
