function [hIntrinsic, mReflex, nlmStiff, VAFT, VAFi, VAFr, tqI, tqR, tqT]= PC_stiffnessID (z, varargin)
% function [hIntrinsic, mReflex, nlmStiff, VAFT]= PC_stiffnessID (z, varargin)
% Parallel cascade stiffness Identificiation 
% Uses: PC_StiffnessID (z,'?') to set optional parmeter name value pairs
% V01-02 REK
% Optional argume name/value pairs are:
% 
% Define options, default values and help 
options={{'reflex_id_method' 'sls' 'method for identification of reflex pathway (MbyZ/TIC/ScanType)'} ...
        {'decimation_ratio' 10 'ratio to decimate data bt'} ...
        {'reflex_irf_len' 1000 'length of reflex IRF in ms'} ...
        {'reflex_nl_ordermax' 8 'maximum order for reflex nonlinearity'} ...
        {'irf_mode' 'pseudo' 'mode for IRF identification {pseudo/full/manual'} ...
        {'id_tolerance' .05 'tolerance on VAF improvement of VAF for main loop'} ...
        {'id_max_iter' 10 'maximum number of iterations for main loop'} ...
 
    };
keyboard
if arg_parse(options,varargin);
    return
end
ts=get(z,'domainIncr');
nlag_reflex_irf= ceil(reflex_irf_len/(ts*1000*decimation_ratio));
mReflex=nlbl;
if strcmpi(reflex_id_method,'hk'),
    set(mReflex,'idMethod','hk','nLagLE',nlag_reflex_irf, 'maxOrderNLE',reflex_nl_ordermax);
elseif strcmpi(reflex_id_method,'sls'),
    set(mReflex,'idMethod','sls', 'displayFlag',0,'nIterMax',10,'nLagLE',nlag_reflex_irf, 'maxOrderNLE',reflex_nl_ordermax);
elseif strcmpi(reflex_id_method,'rect'),
    mReflex=irf;
    set (mReflex,'nSides',1,'nLags',nlag_reflex_irf,'irfIdMethod',irf_mode);
else
    error ('Method not defined');
end
pos = z(:,1); 
vel = ddt(pos);
vel=decimate(vel,decimation_ratio);
vrect=max(vel,0);
z=decimate(z,decimation_ratio); 

nlag_intrinsic=ceil(.04/(ts*decimation_ratio));
hIntrinsic = irf;
hIntrinsic= set(hIntrinsic,'nLags',nlag_intrinsic,'nSides',2);

%set(mReflex,'nLagLE',nlag_reflex_irf, 'maxOrderNLE',reflex_nl_ordermax);

z=z-mean(z);
pos=z(:,1);

tq=z(:,2);
tqR=tq*0;
VAFbest=-1000;
nsamp=length(pos);
for loop=1:id_max_iter,
    %
    % Linear pathway
    %
    tqI=tq-tqR;
    zp=cat(2,pos, tqI);
    %hIntrinsic=irf(hIntrinsic, zp);
    hIntrinsic = nlident(hIntrinsic,zp);
    tqI=nlsim(hIntrinsic,pos);
    %clf;
    %figure(1);
    %plot (hIntrinsic);
    %
    % Nonlinear pahtway with rectifer
    %
    tqR=tq-tqI;
    if strcmpi(reflex_id_method,'rect'),
        zv=cat(2,vrect,tqR);
        mReflex=nlident(mReflex, zv);
        tqR=nlsim(mReflex,vrect);
    else    
    zv=cat(2,vel,tqR);
    mReflex=nlident(mReflex, zv);
    tqR=nlsim(mReflex,vel);
    end
    %figure(2); clf; 
    %plot (mReflex);
    tqT=tqI+tqR;
    tqx=tq(nlag_reflex_irf:nsamp);
    tqTx=tqT(nlag_reflex_irf:nsamp);
    VAFT=vaf(tqx,tqTx);
    VAFT=double(VAFT);
    
    %%HIG 02/23/06
    %VAFi=vaf(tqT(nlag_reflex_irf:nsamp-5),tqI(nlag_reflex_irf:nsamp-5));
    VAFi=double(vaf(tq(nlag_reflex_irf:nsamp-5),tqI(nlag_reflex_irf:nsamp-5)));
    %VAFr=vaf(tqT(nlag_reflex_irf:nsamp-5),tqR(nlag_reflex_irf:nsamp-5));
    VAFr=double(vaf(tq(nlag_reflex_irf:nsamp-5),tqR(nlag_reflex_irf:nsamp-5)));
    
    ss=sprintf('Iteration: %4i \tTotal VAF: %6.3f Delta VAF: %6.3f\n',loop, VAFT,VAFbest-VAFT);
    disp(ss);
    %streamer ([ reflex_id_method ss]);
    %drawnow;
    if (VAFT-VAFbest) < id_tolerance,
        break
    else
        VAFbest=VAFT;
    end
end
%mel=get(mReflex,'elements');
% mddt=nlmsym(sym('ddt(x)'));
nlmStiff=nlm;
set(nlmStiff,'elements', { hIntrinsic ; mReflex } );