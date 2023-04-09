
function date_time = unixt2datestr (unix_time)
% Convert a unix time stamp to a date time str
date_time = datestr(unix_time./86400 + datenum(1970,1,1));
return