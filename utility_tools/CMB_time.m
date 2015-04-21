function output=CMB_time(DBC,option,varargin)
% CMB_TIME Time and date utilities
%  
% Option 
%   current_time : get the date and time based on the db
%                  because the matlab clock function is not 
%                  reliable
%   check_time : verify the time validity 
%   check_date : verify the date validity
% 
% Example 
%   T=CMB_time(DBC,'current_time');
%
%   CMB_time(DBC,'check_time',H24,MI,SS) with H24,MI and SS type=double
% See also 
% 
%% AUTHOR: Yannick Richard 
%% $DATE: 08-Aug-2008 06:59:04 $ 
%% $Revision: 1.4 $ 
%% DEVELOPED: 7.6.0.324 (R2008a) 
%% FILENAME: CMB_time.m 
%Copyright 2002-2005 McGill University This file is part of the CMB Proteomics PipeLine 

switch lower(option)
    %current day time based on the db because datestr(clock,31) is not reliable
    case 'current_time'
        sql_query=['select sysdate(''yyyy-mm-dd HH24:MI:SS'')' ...
                  ' from dual'];
        cursTemp = exec(DBC,sql_query);
        if ~isempty(cursTemp.Message),
            error(cursTemp.Message);
        end
        cf_seq=fetch(cursTemp);
        % we get 2008-12-05 16:52:25.0
        % truncate to have 2008-12-05 16:52:25
        output=cf_seq.Data{1}(1:(length(cf_seq.Data{1})-2));
        close(cf_seq);
        close(cursTemp);

    case 'check_time'
        error(nargchk(3, 5, nargin));
        [hour, minute, second]=varargin{1:3};
        output=istime(hour, minute, second);
        
    case 'check_date'
        error(nargchk(3, 5, nargin));
        [year, month, day]=varargin{1:3};
        output=isdate(year, month, day);
        
    otherwise
        error(['Unknown option: ' option]);
end

end

%%
function t = isdate(time_year, time_month, time_day)
%ISDATE True for valid dates (Gregorian calendar).
%
%   ISDATE(YEAR, MONTH, DAY) returns 1 if input is a valid year-month-date
%   triple and 0 otherwise.  Gregorian calendar is assumed.
%
%   See also ISJDATE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:51:49 +0100
%   E-mail:      pjacklam@online.no
%   URL:  http://home.online.no/~pjacklam/matlab/software/util/timeutil/

   t =   ~imag(time_year) & ~imag(time_month) & ~imag(time_day) ...
       & (time_year == round(time_year)) & (time_month == round(time_month)) ...
       & (time_day == round(time_day)) & (1 <= time_month) & (time_month <= 12) & (1 <= time_day);

   % Since this function might be called at the beginning of other
   % time-related m-files, we make this function independent of all other
   % m-files to avoid infinite recursion.

   days = [ 31 28 31 30 31 30 31 31 30 31 30 31 ];
   is_february = time_month == 2;
   is_leapyear = ( ~rem(time_year, 4) & rem(time_year, 100) ) | ~rem(time_year, 400);
   days_in_month = days(time_month) + ( is_february & is_leapyear );

   t = t & ( time_day <= days_in_month );
end

%%
function t = istime(time_hour,time_minute,time_second)
%ISTIME True for valid times.
%
%   ISTIME(HOUR, MINUTE, SECOND) returns 1 if input is a valid
%   hour-minute-second triple and 0 otherwise.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 12:51:42 +0100
%   E-mail:      pjacklam@online.no
%   URL:   http://home.online.no/~pjacklam/matlab/software/util/timeutil/

   % Hour is an integer, 0 <= hour <= 24.  Minute is an integer, 0 <= minute
   % < 60, with one exception: when hour is 24, minute must be 0.  Second is
   % a real number, 0 <= second < 60, with two exceptions: firstly, when
   % hour is 24, second must be 0; secondly, when hour is 23 and minute is
   % 59, second may be 0 <= second < 61 (to allow for positive leap
   % seconds).  A positive leap second is introduced by letting the last
   % minute of the last hour (of the last day of the month) have 61 seconds.
   % I.e., 23:59:60 <= leap second < 23:59:61 = 00:00:00 the following day.
   
   t =   ~imag(time_hour) & ~imag(time_minute) & ~imag(time_second)     ...  % real
       & time_hour == round(time_hour) & time_minute == round(time_minute)   ...  % integers
       & 0 <= time_hour & 0 <= time_minute & 0 <= time_second           ...  % positive
       & (   ( time_hour <= 23 & time_minute <= 59 & time_second < 60 ) ...  % most times
           | ( time_hour == 23 & time_minute == 59 & time_second < 61 ) ...  % allow leap sec
           | ( time_hour == 24 & time_minute ==  0 & time_second == 0 ) ...  % midnight
         );
end


% Created with NEWFCN.m by Frank González-Morphy  
% Contact...: frank.gonzalez-morphy@mathworks.de  
% ===== EOF ====== [CMB_time.m] ======  
