function y = etvc(x,h,dt,arg1,arg2);
% Ensemble Time Varying Convolution.
% y = etvc(x,h,dt,'typ');
% y = etvc(x,h,dt,n1,n2);
%
% Ensemble Time Varying Convolution.  Ensemble input X is convolved with
% the time varying convolution kernel H,  DT is the intersample interval.  
% X and H have the following structure:
% 
%    X:  ==== Trials ====>      H:  ==== Lag Time ====>
%        ||                         ||
%        ||                         ||
%        ||                         ||
%     Discrete                   Discrete
%       Time                       Time
%        ||                         ||
%        ||                         ||
%        ||                         ||
%        ||                         ||
%        \/                         \/
%
% The filter bounds can be specified using one of two methods, which
% ever is apropriate for the filter, H, being used:
%
% 1. Simple specify 'TYP' where 'TYP' equals 'one' to identify a one sided,
% and 'two' to identify a two sided filter.  The filter length is determined
% automatically from the column dimension of H.
%
% 2. Specify N1, and N2.  Where N1, and N2 are the filter bounds expressed in 
% discrete time. To identify a filter which goes from -2 to 10, use
% N1 = -2, and N2 = 10.  This will yield a 13 point filter.  The filter
% bounds are restricted as follows:  N1 <= 0, and N2 > 0.  The filter
% length (=N2-N1+1) should agree with the column dimension of H. 
%
%
% See also:  TVFIL ETVC_ 

% Copyright 1992-2003, Robert E Kearney and James P Trainor
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

[dummy firl] = size(h);

if nargin == 4
  typ = arg1;
  if strcmp(typ,'two')
    % check if the filter length is odd
    if floor(firl/2) ~= firl/2
      % An odd length two sided filter is split evenly about lag zero.
      % n1 and n2 are in discrete time, they are not indexes
      n1 = -floor(firl/2);
      n2 =  floor(firl/2);
      % the index of the zero lag point in the identified fir will be:
      zero_lag_point = floor(firl/2)+1;
    else   % the filter has even length
      % An even length two sided is split unevenly about lag zero.
      % Positive time gets one more point that negative time.
      % n1 and n2 are in discrete time, they are not indexes.
      n1 = -firl/2 + 1;
      n2 =  firl/2;
    end
  elseif strcmp(typ,'one') 
    % the begining and end of a one sided filter are straight forward,
    % note n1 and n2 are discrete times, they are not indices
    n1 = 0;
    n2 = firl-1;
  else
    error('TYP can only be ''one'' or ''two''')
  end
elseif nargin == 5
  n1 = arg1;
  n2 = arg2;
  if ~(n1 <= 0  &  n2 > 0  &  (n2-n1+1) == firl)
    error('improper filter bounds, n1,n2');
  end
else
  error('invalid number of input arguments')
end


% now call the mex file which does the actual convolution
y = etvc_mex(x,h,dt,n1,n2);

% end of file etvc.m


