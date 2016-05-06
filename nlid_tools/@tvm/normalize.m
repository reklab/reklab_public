function [Hnorm] = normalize(H, mode)

% [Hnorm] = normalize(H) 
%
% Mode = 1: Changes the gain of the static nonlinearity and of the dynamic linear system 
% at every time i such that the overall gain remains the same but the dc gain 
% of the linear system is 1.  
%
% Mode = 2: Changes the gain of the static nonlinearity and of the dynamic linear system 
% at every time i such that the overall gain remains the same but the dc gain 
% of the polynomial is 1.  

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

% Check that number of input arguments is valid.
if nargin ==1,
    mode=1;
end
num_args = nargin;


MT=get(H,'Model_Type');
if ~strcmp(MT,'nlbl'),
    error('model type not support by normalize');
end
Hnorm=H;
hd=get(H,'data');
for i=1:length(hd),
    mtemp=hd{i};
    hel=get(mtemp,'elements');
    p=hel{1};
    I=hel{2};
    dt=get(I,'domainincr');
    Id=get(I,'data');
    Pc=get(p,'coef');
    if mode==1,
        
        factor=sum(Id)*dt;
        
    elseif mode==2,
            factor =1/Pc(2);
        end
        Idnew=Id/factor;
        PcNew=Pc*factor;
        set(p,'coef',PcNew);
        set(I,'data',Idnew);
        hel{1}=p;
        hel{2}=I;
        set(mtemp,'elements',hel);
        hnew{i,1}=mtemp;
    end
    set(Hnorm,'Data',hnew);
    
    
    
    
    %
    
    % end @tv_hammerstein\normalize.m
    
    
