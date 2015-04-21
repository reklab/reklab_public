function struct_table (fname,S,TableTitle, cf)
%struct_table - output a structure to a html table
% $Revision: 1.1.1.1 $
% Usage:
%   struct_disp (fname, S, TableTitle)
%  fname - output filename
%  S - structure to display
%   cw - vector of colmn widths or scalar specifying all widths
% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
if length(fname)==0,
    fid=1;
else
    fid=fopen(fname,'w');
end
names = fieldnames(S);
%
if iscell(TableTitle),
    for i=1:length(TableTitle),
        fprintf (fid,'<center><h1> %s </h1></center><p>',TableTitle{i});
    end
else
    fprintf (fid,'<center><h1> %s </h1></center><p>',TableTitle);
end
fprintf (fid,'<table>');
fprintf(fid,'<tr>\n');
for i=1:length(names),
    fprintf(fid,'<th> %s </th>',names{i}); 
end
fprintf(fid,'</tr>');

for i=1:length(S),
    fprintf(fid,'<tr>');
    for j=1:length(names),
        tblcol(fid,S(i).(names{j}));
    end
    fprintf(fid,'</tr>');
end
fprintf(fid,'</table>');
if fid > 2,
    fclose(fid);
end

function tblcol(fid, el);
if ischar(el),
    fprintf(fid,'<td> %s </td>\n', el); 
elseif isnumeric(el),
    fprintf(fid,'<td> %s <//td>\n', num2str(el)); 
elseif islogical(el),
    if el,
        fprintf(fid,'<td> true <//td>\n'); 
    else
        fprintf(fid,'<td> false <//td>\n'); 
    end
elseif ischar(el{1}),
    fprintf(fid,'<td> %s </td>\n', el{1}); 
    
end




