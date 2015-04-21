function [R, V, yp] = nlid_resid( M, z, varargin);
% NLID_RESID - compute and display prediction error in model output.
% Usage:
%   [R, V, yp] = nlid_resid( M, z, plotFlag);
%       R - residuals
%       V - variance accoutned for
%      yp - pedicted otuput
%
% M - model
% z - innout output dats


% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt
options={{'plotflag' true 'Plot results '} ...
         {'choplen' 0 'length of transiet to ignore at begining and end of response'} ...
      
     };
 if arg_parse(options,varargin);
     return
 end
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
    plot (cat(2,y,yp),'plotmode','Super');
    title('Superimposed');
    subplot (4,1,4);
    plot (R);
    T= ('Residuals' );
    title(T);
end
  disp(Vt);

