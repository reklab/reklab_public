function [token]=html_token ( S, start_string, stop_string, string_flag, ntoken)
% [token]=html_token ( S, start_string, stop_string, string_flag, ntoken)
% Returns tokens  in string delimited by start_string and strop string
% Default is to return as a cell array
% if string flag ==1 returns a string
% $Revisoins: $
%2002 10 15 REKhankde zero length strings
% 2002 12 18 REK Modify to handle stop string of length > 1
if strcmp(start_string, stop_string); 
    error('html_token: start and stop strings must be distinct');
end
if nargin <4,
    string_flag=0;
end
if nargin <5,
    ntoken=0;
end
istart=strfind(S, start_string);
istop=strfind(S,stop_string);
token={};
j=1;
for i=1:length(istart),
    itstart= istart(i)+length(start_string); 
    itstop= min(istop(istop>=itstart)) - 1;
    if itstop==0,
        break
    elseif itstop<itstart,
        token{j}='';
        if ntoken==j,
            break
        end
          j=j+1;    
    else
        token{j}=S(itstart:itstop);
        if ntoken==j,
            break
        end
          j=j+1;
    end
  
end

if string_flag,
    token=char(token);
end

