function nlout = hammerhk ( z, nl, p );
% Identification of Hammerstein cascade using the Hunter-Korenberg
% iterative cross-correlation method.
%
% See: I.W. Hunter and M.J. Korenberg, The identification of nonlinear
%      biological systems: Wiener and Hammerstein cascade models.
%      Biological Cybernetics, 55:135-144, 1986.

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../../copying.txt and ../../gpl.txt

%
% 23 May 2002 - Normalization needs work and will probably not work
% for HiPass systems

assign(p) % Assign paramter values in workspace
u=z(:,1);
y=z(:,2);
yd=double(y);
xEst=u;
outputSS=yd'*yd; % sum of squares of output
oldSS = -inf;
newSS = 1;
yp=y;
invi=irf;
set(invi,'nSides',2,'nLags',nLagLE,'irfPseudoInvMode','auto');
i=irf;
set(i,'nSides',1,'nLags',nLagLE,'irfPseudoInvMode','auto');
flag=1;
h_final=nan;
while flag,
    oldSS=newSS;
    z1 = cat (2,y,xEst);
    z1=z1-mean(z1);
    inv_h = nlident(invi, z1);
    halflen = floor(length(inv_h)/2);
    xEst = nlsim (inv_h, y);
    z2=extract(cat(2,u,xEst),halflen);
    t=nl.elements{1};
    m_est = polynom(t,z2, 'polyType','tcheb','polyOrderMax',maxOrderNLE );
    xEst = nlsim (m_est, u );
    z3= cat (2, xEst,y );
    h_est = nlident(i, z3);
    %% Adjust DC terms for low pass systems
    N=nlbl('elements',{m_est h_est});
    N=normGainLE(N);
    e=N.elements;
    m_est=e{1}; h_est=e{2};
    pc=m_est.polyCoef;
    pc(1)=0;set(m_est,'polyCoef',pc);
    N=nlbl('elements',{m_est h_est});
    yp=nlsim(N,u);
    delt = mean(y)-mean(yp);
    pc(1)=double(delt);set(m_est,'polyCoef',pc);
    N=nlbl('elements',{m_est h_est});
    yp=nlsim(N,u);
    
    %%
    
    currentError = double(yp-y); 
    newSS=currentError'*currentError/outputSS;  % sub of square of current error
    deltSS=oldSS-newSS;
    if (displayFlag)
        fprintf ('Old SSE %6.4f New SSE %6.4f Delta SSE %6.4f \n', oldSS, newSS,  deltSS);
    end
    if (deltSS)> threshNSE | isnan (h_final)
        h_final=h_est;
        m_final=m_est;
    else
        break
    end
end
%
% Adjustment for DC offset
%
h_est=h_final;
set(h_est,'comment','Linear Element');
m_est=m_final;
set(m_est,'comment','Static Nonlinearity');
elements = { m_est h_est };
Ts=z.domainIncr;
set (nl,'elements',elements);
nlout=nl;

