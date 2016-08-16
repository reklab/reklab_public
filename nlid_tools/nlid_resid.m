function [R, V, yp] = nlid_resid( M, z, varargin);
% NLID_RESID - compute and display prediction error in model output.
% Usage:
%   [R, V, yp] = nlid_resid( M, z, plotFlag);
%       R - residuals
%       V - variance accoutned for
%      yp - pedicted otuput
%
% M - model
% z - inpout output dats
% Multiple realizations are handled as follows
% NM - number of model realizationa
% NZ - number of data realizations
%  If NM==1 & NZ > 1
%     errors is computed for each realization.
%  IF NM > 1 7 NZ > 1
%     error is computer for each model wioth the data
% If NM > 1 NZ > 1 and NM=NZ
%  Run ech model for the associated data set
% if NM>1 & Nz>1 & NM ~= NZ
%  error
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU

[ nSampM, nChanM, nRealM]=size(M);
[ nSampZ, nChanZ, nRealZ]=size(z);

if (nRealM>1) & (nRealZ>1) & (nRealM ~= nRealZ)
    error (' There is a mismatch in the number of realizations for the model and data');
end
options={{'plotflag' true 'Plot results '} ...
    {'choplen' 0 'length of transiet to ignore at begining and end of response'} ...
    };
if arg_parse(options,varargin);
    return
end
iReal=1;

for iRealM=1:nRealM,
    for iRealZ=1:nRealZ,      
        if iRealM==1 & iRealZ==1,
            [R, V, yp] = nlid_resid_sub ( M, z, plotflag, choplen);
        elseif nRealM==1 && iRealZ>1,
            iReal=iReal+1;
            [R(:,:,iReal), V(iReal), yp(:,:,iReal) ] = nlid_resid_sub ( M, z(:,:,iRealZ), plotflag, choplen);
        elseif iRealM>1 && nRealZ==1,
            iReal=iReal+1;
            [R(:,:,iReal), V(iReal), yp(:,:,iReal) ]=  nlid_resid_sub ( M(:,:,iRealM), z, plotflag, choplen);
        elseif iRealM > 1 && iRealZ > 1 && iRealM==iRealZ,
            iReal=iReal+1;
            [R(:,:,iReal), V(iReal), yp(:,:,iReal) ]= nlid_resid_sub ( M(:,:,iRealM), z(:,:,iRealZ), plotflag, choplen);
        end
    end
end

end

function [R, V, yp] = nlid_resid_sub ( M, z, plotflag, choplen);

if isa(M,'polynom');
    nin=M.nInputs;
    x=z(:,1:nin);
    y=z(:,nin+1);
else
    x=z(:,1);
    y=z(:,2);
    
end

yp= nlsim(M,x); yp=yp(:,1);
% Get rid of transients
y=chop(y,choplen);
yp=chop(yp,choplen);

R=y-yp;

comment=M.comment;
set(R,'comment',['Residuals of '  comment]);
V=vaf(y,yp);
V=double(V);
Vt=['%VAF = ' num2str(chop(double(V),4))];

if plotflag,
    
    subplot (4,1,1);
    plot(y);
    title('Observed');
    subplot (4,1,2);
    plot (yp);
    
    
    title(['Predicted ' Vt]);
    
    subplot (4,1,3);
    plot (y);
    h=line(yp); set(h,'color','r'); 
    title('Superimposed');
    subplot (4,1,4);
    plot (R);
    T= ('Residuals' );
    title(T);
end
disp(Vt);

end