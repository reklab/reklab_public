function params = lm_sls_parameters(tp)
% function to set up training parameters for Levenberg-Marquardt iterations.
%
%        tp is a vector containing training parameters
% 
% tp = [max_its, threshold, accel, decel, delta, display];
% if some entries are not specified, or passed as nans, the following 
%   defaults are assigned.
%
%   max_its: 20 (iterations)
%   threshold: 0 (normalized mean squared error).
%   accel:  0.8  (ridge is multiplied by 0.8 following a successful update).
%   decel:  2    (ridge is multiplied by 2 following a failed update).
%   delta:  10   (initial size of the ridge added to the Hessian)
%   display: 1   (binary, display intermediate results)

% Copyright 2000-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

tp_defaults = [20 0 0.8 2 10 0]'; 

tp = tp(:);
params = nan*zeros(6,1);

params(1:length(tp)) = tp;
inds = isnan(params);
params(inds) = tp_defaults(inds);








