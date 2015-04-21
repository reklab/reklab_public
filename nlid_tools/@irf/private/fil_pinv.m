function [irf,bound] = fil_pinv(u,y,hlen,sides,level,mode);
% nonparamtric filter identification
% solves for the nonparametric filter between a given
% input and output using correllation techniques, and 
% Toeplitz matrix (pseudo) inversion.
%
% Syntax: [nfil,bounds,mdl]= fil2(u,y,numlags,numsides,level,mode)
%
%
%    u        : input to the unknown system
%    y        : self-explanatory I hope
%    numlags  : the number of lags to calculate the linear filter.  Note,
%               for 2 sided filters, this will be rouded up to the nearest
%               ODD number, so that the filter is symmetric about zero lag.
%    numsides:  determine a causal (1 sided, default) or noncausal 
%		(2 sides) response.  The noncausal response is
%		computed by shifting the input forward in time
%		by numlags/2 points, hence, for a two sided
%		identification, the following condition must
%		hold	numpts < size(x) - numlags/2
%    level    : if confidence bounds are desired, this is the desired
%               confidence level (1 - 99%, default 95%)
%    mode     : Method used to determine the rank of the identifiable
%               subspace used in the low-rank projection 
%               'full', 'auto', 'manual'  (default 'auto')
%               if 0<mode<numlags+1, then mode is used as the pseudo-inverse
%               order.
%
%
%    nfil     : the identified nonparametric filter.
%    bounds   : confidence bounds on the impulse response
%    mdl      : Minimum Descriptor Length cost function used to determine
%               the rank of the identifiable subspace.  (Returns zeros if
%               mode is 'full')
%

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

global dtw_pinv_input dtw_pinv_output dtw_pinv_numlags dtw_pinv_numsides
global dtw_pinv_order dtw_pinv_level

global dtw_pinv_irf dtw_pinv_bound
%global dtw_pinv_IRFs dtw_pinv_Bounds 
global dtw_pinv_done

global dtw_pinv_done
global dtw_mdls

%
% Parse input parameters
%

if nargin < 6
  mode = 'auto';
  if nargin < 5 
    level = 95;
    if nargin < 4
      sides = 1;
    end
  end
end

if length (sides)==0,
  sides=1
end
if length(level)==0,
  level=95;
end

dtw_pinv_input = u;
dtw_pinv_output = y;
dtw_pinv_numlags = hlen;
dtw_pinv_numsides = sides;
dtw_pinv_level = level;
compute_pinv('init');



if mode(1) == 'a'
  [val, dtw_pinv_order] = min(dtw_mdls);
  compute_pinv('irf');
elseif mode(1) == 'f'
  dtw_pinv_order = hlen;
  compute_pinv('irf');
elseif (~isstr(mode) & mode > 0 & mode < hlen+1)
  dtw_pinv_order = mode;
  compute_pinv('irf');
elseif mode(1) == 'm'
  [val, dtw_pinv_order] = min(dtw_mdls);
  pinv_gui('init');
  %  use a tight loop, and wait for choose_irf to complete.
  done = 0;  
  while done == 0;
    done = dtw_pinv_done;
    pause(0.25);
  end
else
  error('unsupported mode');
end  


irf = dtw_pinv_irf;
bound = dtw_pinv_bound;
