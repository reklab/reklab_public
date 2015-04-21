function struct2txt(S, fileName, varargin)
%struct_table - output a structure to a csv file 
% $Revision: 1.1 $
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
    {'title' 'default'} ...
    {'browser','-browser'} ...
    {'delimiter', ','} ...
    {'display_flag' true} ...
    {'display_format' 'column'} ...
     {'subtitle' {} } ...
     {'separatorstring' '<br>' } ...
    {'permission' 'wt' }};
arg_parse(options, varargin);
 if ~iscell(S),
     temp={S};
     S=temp;
 end
 
fid = openOutputFile (fileName, permission) ;
for iTable=1:length(S),
    writeHeader (S{iTable}, fid, delimiter); 
    writeData (S{iTable},fid, delimiter);
end
fclose (fid) 

disp(fid)

end

function fid = openOutputFile (fileName, permission)
fid=fopen(fileName,permission);
if fid <1,
    error('Error opening file');
end
end

function writeHeader (S, fid, delimiter) 
% Write header line
headers=fieldnames(S);
for i = 1:length(headers)
   fprintf(fid,headers{i});
    if i < length(headers)
        fprintf(fid,delimiter);
    end
end
fprintf(fid,'\n');
end

function writeData (S, fid, delimiter)
% Determine the data types for each field
header=fieldnames(S);
nVal=length(S.(header{1}));
nField=length(header);
for iField = 1:length(header)
    if isnumeric(S.(header{iField})(1))
        cellFormat{iField} = 'numeric';
    elseif iscell(S.(header{iField})(1))
        cellFormat{iField} = 'cell';
    elseif ischar(S.(header{iField})(1))
         cellFormat{iField} = 'char';
    end
end

% Write the data

for iVal = 1:nVal
    for j = 1:nField
        switch cellFormat{j}
            case 'numeric'
                fprintf(fid,'%s',num2str(S.(header{j})(iVal)));
            case 'cell'
                fprintf(fid,'%s',S.(header{j}){iVal});
        end

        if j < length(header)
            fprintf(fid,delimiter);
        end
    end
    fprintf(fid,'\n');
end
end
