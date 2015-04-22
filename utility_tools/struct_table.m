function struct_table (fname,S,TableTitle, varargin)
%struct_table - output a structure to a html table
% $Revision: 1.26 $
% Usage:
%   struct_disp (fname, S, TableTitle)
%  fname - output filename
%  S - structure to display. S may be a structure or cel array of
%    strucutures in which case each element in the cell array is displayed as
%    a separate subtable. 
%  varargin - options defined as name/value pairs
%     browser - browser option
%           - specifies matlab options to Web 
% table_format - normal/transpose
%            normal - structures of format S(n).name
%            transpose - structres of format S.name(n)
%     display_format [column/row]
%                    column - dispay independnet variables in columns
%                    row - display independnet variable in rows
%     display_flag - pop up browser with results
%     permission - permission for file create
%                     'w' - write
%                     'a' - append
%     seperatorstring - string to separate multiple entires in a cell.
%                   defauilt is <br> whjich will generate a newline
%                   use ',' to generate comma separated lists
%     subtitle - cell array of talbe subtitles. One element for each table
%      
% struct_disp(S)
%      displays S in web browser.
% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt
%%
options = { {'table_format' 'normal'}  ...
    {'browser','-browser'} ...
    {'display_flag' true} ...
    {'display_format' 'column'} ...
     {'subtitle' {} } ...
     { 'separatorstring' '<br>' } ...
    {'permission' 'w' }};
arg_parse(options, varargin);
if nargin==1,
    S=fname;
    fname=fullfile(pwd,'temp.html');
    TableTitle='Default';
end
if length(fname)==0,
   fname=tempname;
end
    fid=fopen(fname,permission);
    if fid <0,
        error(['Error openining file:' fname]);
    end
%
if iscell(TableTitle),
    for i=1:length(TableTitle),
        fprintf (fid,'<center><h1> %s </h1></center><p>',TableTitle{i});
    end
else
    fprintf (fid,'<center><h1> %s </h1></center><p>',TableTitle);
end

%% Display table variables by column
if ~iscell(S),
    S={S};
end
for iCell=1:length(S),
    names = fieldnames(S{iCell});
    if ~isempty(subtitle),
        fprintf (fid,'<h1> %s </h1><p>',subtitle{iCell});
    end
    fprintf (fid,'<table rules=''all''>');
    fprintf(fid,'<tr>\n');
    switch display_format
        case 'column'
            for i=1:length(names),
                fprintf(fid,'<th> %s </th>',names{i});
            end
            fprintf(fid,'</tr>');
            switch table_format
                %% format=transpose for structures of form S(i).element;
                case 'normal'
                    for i=1:length(S{iCell}),
                        fprintf(fid,'<tr>');
                        for j=1:length(names),
                            tblcol(fid,S{iCell}(i).(names{j}),separatorstring);
                        end
                        fprintf(fid,'</tr>');

                    end

                    %% format=transpose for structures of fomr S.element(i);
                case 'transpose'
                    for i=1:length(S{iCell}),
                        fprintf(fid,'<tr>');
                        for k=1:length (S{iCell}(i).(names{1})),
                            for j=1:length(names),
                                tblcol(fid,S{iCell}(i).(names{j})(k), separatorstring);
                            end
                            fprintf(fid,'</tr>');
                        end
                    end
                otherwise
                    error('Bad value for table_format');
            end
        case 'row'
            nField=length(names);
            for iField=1:nField
                fprintf(fid,'<tr><td> %s </td>',names{iField});
                tblelement(fid,S{iCell}.(names{iField}));
                fprintf(fid,'</tr>');
            end

        otherwise
            error (['Bad value for display_format:' display_format])
    end
      fprintf(fid,'<tr></tr>');
      fprintf(fid,'</table>');


end
if fid > 2,
    fclose(fid);
end


if (display_flag),
    web(fname, browser);
end

%%


function tblcol(fid, el, separatorString);
%% Print a column of a table
if iscell(el) & length(el)==1,
    el=el{1};
end
if ischar(el),
    fprintf(fid,'<td valign=top> %s </td>\n', el);
else
    fprintf(fid,'<td valign=top>');
    for i=1:length(el),
        if i>1,
              fprintf(fid,'%s' , separatorString);
        end
        if isnumeric(el),
            fprintf(fid,'%s' , num2str(el(i)));
        elseif islogical(el),
            if el(i),
                fprintf(fid,'true' );
            else
                fprintf(fid,'false ' );
            end
        elseif iscell(el),
            if ischar(el{i}),
                fprintf(fid,'%s'  , el{i});
            else
                fprintf(fid,'%s' , num2str(el{i}));
            end
        end

           
    end
    fprintf(fid,' </td>\n');


end


%% Print one element of a column

function tblelement(fid, el);
%% Print an element of a table
if ischar(el),
    fprintf(fid,'<td valign=top> %s </td>\n', el);
else
    fprintf(fid,'<td valign=top>');
    for i=1:length(el),
        if isnumeric(el),
            fprintf(fid,' %s ;', num2str(el(i)));
        elseif islogical(el),
            if el(i),
                fprintf(fid,' true ;');
            else
                fprintf(fid,' false ; ');
            end
        elseif iscell(el),
            if ischar(el{i}),
                fprintf(fid,' %s ; ', el{i});
            else
                fprintf(fid,' %s ;', num2str(el{i}));
            end
        end
    end
    fprintf(fid,' </td>\n');
end

