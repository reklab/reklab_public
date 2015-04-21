function struct_print (file,S)
% print table with heading
% Usage:
%
%   table_print (file,table,names,cw_in,form_in,prec_in)
%
%       file - output file (dump to screen if zero length)
%	table 
%	names 
%	cw_in integer vector - column width
%	form_in  - colum format ( integer(i), string(s),float(f))
%       prec_iin - column precision
%
%

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

names = fieldnames (S);
ncol=length(names);
nrow=length(S);

for i=1:80,
    blank_line(i)=' ';
end
for i=1:ncol,
    cw(i)=10;
end

header=blank_line;
current_col =0;
[m,n]=size(names);
for i=1:m,
    name=deblank(names{i});
    ilast=current_col+cw(i);
    n1 = min(length(name),cw(i));
    ifirst=ilast-n1+1;
    header(ifirst:ilast)=name(1:n1);
    current_col=ilast+1;
end
if(length(file)==0),
    disp(header);
else,
    fprintf(file,header);
    fprintf(file,'\n');
end



for jrow=1:nrow,
    current_col=0;
    line=blank_line;
    for icol=1:ncol,
        x = getfield (S, {jrow},names{icol});
        if(isnumeric(x)),
            s=num2str(x);
        else	
            s=sprintf(setstr(x));
        end
        n1=length(s);
        if (n1 >cw(i)),
            s='*';
            n1=1;
        end
        ilast = current_col+cw(i);
        ifirst = ilast-n1+1;
        line(ifirst:ilast)=s(1:n1);
        current_col=ilast+1;
    end
    if (length(file)==0),
        disp(line);
    else,
        fprintf (file,line);
        fprintf (file,'\n');
    end
end







