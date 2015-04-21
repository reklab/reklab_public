function Q=struct_disp (S, cw, varargin)
%struct_disp - display  struct with heading $Revision $
% $revision: $
% Usage:
%   struct_print (S, cw)
%  S - structre to display
% cw - vector of colmn widths or scalar specifying all widths
% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

options = {{ 'file' 'unknown'}};
arg_parse(options, varargin);

if strcmp(file,'unknown')
    file=1;  % display on screen
end
if iscell(S) & length(S)==1,
    newS=S{1};
    S=newS;
end
if length(S)==0 | isempty(struct2cell(S)),
    disp('empty');
    return
end

names = fieldnames (S);
ncol=length(names);
nrow=length(S);


if nargin==1,
    for i=1:ncol,
        cw(i)=10;
    end
end

if length(cw)==1,
    colwidth=cw;
    for i=1:ncol
        cw(i)=colwidth;
    end
end


for i=1:232,
    blank_line(i)=' ';
end

header=blank_line;
current_col =1;
[m,n]=size(names);
for i=1:m,
    name=deblank(names{i});
    n1 = min(length(name),cw(i));
    ifirst=current_col;
    ilast=current_col+n1-1;
    header(ifirst:ilast)=name(1:n1);
    current_col=ifirst+cw(i)+1;
end
if(length(file)==0),
    disp(header);
else,
    fprintf(file,header);
    fprintf(file,'\n');
end

for jrow=1:nrow,
    current_col=1;
    line=blank_line;
    for icol=1:ncol,
        x = S(jrow).(names{icol});
        if(isnumeric(x)),
            s=num2str(x);
        elseif islogical(x),
            s=num2str(x);
        elseif iscell(x);
            s=' ';
        else
            s=sprintf(x);
        end
        n1=length(s);
        if (n1 >cw(icol)),
            s=s(1:cw(icol));
            n1=cw(icol);
        end
        ifirst=current_col;
        ilast=current_col+n1-1;
        line(ifirst:ilast)=s(1:n1);
        current_col=ifirst+cw(icol)+1;
    end
    if (length(file)==0),
        disp(line);
    else,
        fprintf (file,line);
        fprintf (file,'\n');
    end
end