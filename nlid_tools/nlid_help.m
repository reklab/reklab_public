function nlid_help(option)
% NLID_HELP - generate help on the classes for nlid_tools
% option
%  '  buildHelp - build extended html help file' ...
%    classList - list of all classes for nlid_tools' ...
%     help = display detailed help file' ...
%     classname - help on methods for class' ... 
%  
% Usage
%  nlid_help(option)
% 
% Example
% list all calss for nlid_tools
%   nlid_help ('classList');
% 
% display extend web hep on nlid_tools
%   nlid_help('help');
%
% help on methods for nldat 
%
%   nlid_help ('nldat'); 
%
% See also
%
%% AUTHOR: Robert Kearney
%% $DATE: 22-Mar-2007 11:45:06 $
%% $Revision: 1.2 $
%% DEVELOPED: 7.3.0.267 (R2006b)
%% FILENAME: nlid_help.m
%Copyright 2002-2005 McGill University This file is part of the CMB
%Proteomics PipeLine
helpFileName='nlid_tools.html'; 
optionList = { 'Options are:' ...
    '  buildHelp - build extended html help file' ...
    '  classList - list of all classes for nlid_tools' ...
    '  help = display detailed help file' ...
    '  classname - help on methods for class' ...
    };
if nargin==0,
    for i=1:length(optionList),
        disp(optionList{i})
    end
    return
end
nlidDir=fileparts(which('nlid_help'));
switch lower(option)
    case 'classlist'
        classListHelp(nlidDir);
    case 'help'
        w=which(helpFileName);
        if isempty(w),
            disp('Help file does not exist. Building');
            nlid_help('buildHelp');
        end
        web(helpFileName); 
    case 'buildhelp'
        help2doc(nlidDir)
    otherwise
        classHelp(nlidDir, option);
end
end

function classListHelp(nlidDir);
% Classess list classes
disp(['nlid_tools. Defined classes']);
D=dir(nlidDir);
iDir=find([D.isdir]);
classD=D(iDir);
for i=1:length(classD),
    className=classD(i).name;
    if strcmp(className(1),'@'),
        tmpHelp=help ([className '/' className(2:end)]);
        disp(['  ' className ':' strtok(tmpHelp,char(10))])
    end
end
end

function classHelp( nlidDir, className);
% Class list classes
disp(['nlid_tools. Methods for class: ' className]);
D=dir([nlidDir '/@' className]);
iFile=find(~[D.isdir]);
classD=D(iFile);
for i=1:length(classD),
    fileName=classD(i).name;
    [a,b,c]=fileparts(fileName);
    tmpHelp=help ([className '/' b]);
    disp([ '  ' className '/' b ':' strtok(tmpHelp,char(10))])
end
end

function [index, count]=help2doc(targetdir, fh1, fh2, index, count, h)
% HELP2DOC copies help text for all m-files in a folder (and subfolders)
% to an HTML file
% 
% Examples:
% HELP2DOC()
% HELP2DOC(foldername)
% 
% If foldername is not specified, the user will be prompted for it.
% 
% HELP2DOC copies the help text from all m-files in foldername together
% with any in subfolders. The subfolders are searched to any depth 
% recursively. 
% Text for each m-file is bookmarked internally in the output file.
% HELP2DOC also adds a contents list and an alphabetical index to the start
% of the file. These are hyperlinked internally to the help text.
% 
%
% Ackowledgements: HELP2DOC contains a modified version of the MATLAB
%                  R2006b helpwin function
%
% Revisions: makehelphyper included from matlab helptools
%
% See also helpwin
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 02/07
%         King's College London 2006
% -------------------------------------------------------------------------
%


% If none supplied, ask for the target directory
if nargin==0
    targetdir=uigetdir();
end

% If invoked from the command line, need to open the HTML file for output..
if nargin<2
     % This file receives the contents, and finally the index and contents of
    % help2doc.tmp
    % [filename pathname]=uiputfile('*.html', 'Select file for output', [targetdir '.html']);
    filename='nlid_tools.html';
    pathname=fileparts(which('nlid_help'));
    fh2=fopen(fullfile(pathname,filename),'w+');
    % Temporary file for help text
    fh1=fopen('help2doc.tmp','w+');
    if fh1<0 || fh2<0
        error('help2doc: Failed to open file for output');
    end
    fprintf(fh2, '<html><title>%s generated from help2doc</title>', filename);
    % Contents bookmark & title
    fprintf(fh2, '<p><a href="#Index">[Go To Alphabetical Index]</a></p>');
    fprintf(fh2, '<a name=Contents></a>');
    title = '<p><table width="100%" border=0 cellspacing=0 bgcolor=d0d0f0><tr>';
    title = [title '<td align=center>&nbsp;<b>Contents</b></td></tr></table>'];
    fprintf(fh2, '%s<p>', title);
    h=waitbar(0,'help2doc');
    set(h,'Name','help2doc')
%...and initialize the index
    index=cell(99999,1);
    count=1;
end

% Process the m-files
d=dir(targetdir);
d=dirsort(d);

if length(dir([targetdir filesep '*.m']))>0
    title = '<table width="100%" border=0 cellspacing=0 bgcolor=d0d0f0><tr>';
    title = [title '<td>&nbsp;<b>' targetdir '</b></td></tr></table>'];
    fprintf(fh2, '%s', title);
end

for i=1:length(d)
    [pathname filename extension]=fileparts(d(i).name);
    % Search only m-files and ignore contents.m files
    if strcmpi(extension,'.m')==1 && strcmpi('contents',filename)==0
        waitbar(i/length(d), h, ['Processing file....' filename]);
        file=[targetdir filesep filename extension];
        out=helpwin2(file, count);
        if ~isempty(out)
            index{count}=filename;
            % File 1 - temporarily receives help text
            fprintf(fh1, '<p><a name=b%s></a> %s</p></b>', num2str(count),out);
            fprintf(fh1, '<p><a href="#Contents">[Contents]</a><a href="#Index">[Alphabetical Index]</a></p>');
            % File 2 - receives contents list
            fprintf(fh2, '<p><a href="#b%s">[%5s]%75s</a></p>',...
                num2str(count), num2str(count), filename);
            count=count+1;
        end
    end
end

% Call help2doc recursively for each subfolder
for i=1:length(d)
    if d(i).isdir==true
        folder=[targetdir filesep d(i).name];
        [index, count]=help2doc(folder, fh1, fh2, index, count, h);
    end
end

% If this instance was called by the user, finish off by closing the HTML
% file and deleting the message box
if nargin<2
    % Write the index of functions
    fprintf(fh2, '<html><code><pre><body>');
    % Bookmark the index
    fprintf(fh2, '<a name=Index></a>');
    % Put a title banner up
    title = '<table width="100%" border=0 cellspacing=0 bgcolor=d0d0f0><tr>';
    title = [title '<td align=center>&nbsp;<b>Alphabetical Index of Functions</b></td></tr></table>'];
    fprintf(fh2, '%s', title);
    % Sort the index alphabetically
    index=index(1:count-1);
    [str idx]=sort(lower(index));
    index=index(idx);
    % Genrate the list with hyperlinks
    for i=1:length(index)
        fprintf(fh2, '<p><a href="#b%s">%-70s</a>%10s</a></p>', ...
            num2str(idx(i)), index{i}, num2str(idx(i)));
    end
    % Copy the help text and tidy up
    fseek(fh1, 0, 'bof');
    title = '<p><table width="100%" border=0 cellspacing=0 bgcolor=d0d0f0><tr>';
    title = [title '<td align=center>&nbsp;<b>Function Help</b></td></tr></table>'];
    fprintf(fh2, '%s<p>', title);
    buffer=fread(fh1, Inf, 'uint8');
    fwrite(fh2, buffer, 'uint8');
    fprintf(fh2,'%s',sprintf('</body></pre></code></html>'));
    name=fopen(fh1);
    fclose(fh1);
    delete(name);
    
    name=fopen(fh2);
    fclose(fh2);
    delete(h);
    web(name, '-helpbrowser');
end


return
end


%--------------------------------------------------------------------------
function d=dirsort(d)
%--------------------------------------------------------------------------
% Sorts the  lists alphabetically by file/folder name.
% Note that names are cast to lower case before sorting
name=cell(1,length(d)-2);
d=d(3:end);
j=1;
for i=1:length(d)
    name{j}=lower(d(i).name);
    j=j+1;
end
[str idx]=sort(name);
d=d(idx);
return
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function outStr=helpwin2(topic, count)
%--------------------------------------------------------------------------
% HELPWIN2 is a modified copy of MATLAB's HELPWIN used by help2doc
%
% HELPWIN Online help displayed in the Help window
%   HELPWIN TOPIC displays the help text for the specified TOPIC inside the
%   Help window.  Links are created to functions referenced in the 'See Also'
%   line of the help text.
%
%   HELPWIN(HELP_STR,TITLE) displays the string HELP_STR in the help
%   window.  HELP_STR may be passed in as a string with each line separated
%   by carriage returns, a column vector cell array of strings with each cell
%   (row) representing a line or as a string matrix with each row representing
%   a line.  The optional string TITLE will appear in the title banner.
%
%   HELPWIN({TITLE1 HELP_STR1;TITLE2 HELP_STR2;...},PAGE) displays one page
%   of multi-page help text.  Note: this calling sequence is deprecated and
%   is provided only for compatibility with previous versions of HELPWIN.
%   The multi-page help text is passed in as a
%   cell array of strings or cells containing TITLE and HELP_STR pairs.
%   Each row of the multi-page help text cell array (dimensioned number of
%   pages by 2) consists of a title string paired with a string, cell array
%   or string matrix of help text.  The second argument PAGE is a string
%   which must match one of the TITLE entries in the multi-page help text.
%   The matching TITLE represents the page that is to be displayed first.
%   If no second argument is given, the first page is displayed.
%
%   HELPWIN displays the default topic list in the Help browser.
%
%   See also HELP, DOC, LOOKFOR, WHAT, WHICH, DIR, MORE.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.2 $ $Date: 2007-03-28 15:22:51 $


[pathname fcnName]=fileparts(topic);
helpStr = help(topic);

% Handle characters that are special to HTML
helpStr = fixsymbols(helpStr);

% Is there help for this topic?
if isempty(helpStr)
    outStr=[];
    return
end

% Highlight occurrences of the function name (unless it has ".", "/", or "\"
% immediately after it).
if ~isempty(fcnName) && ~strcmpi(fcnName,'matlab')
    upperFcnName = upper(fcnName);
    lowerFcnName = lower(fcnName);
    if strcmp(lowerFcnName,fcnName) == 0
        helpStr = regexprep(helpStr, ['\<' fcnName '\>(?![\.\/\\])'], ['<b>' fcnName '</b>']);
    else
        helpStr = regexprep(helpStr, ['\<' upperFcnName '\>(?![\.\/\\])'], ['<b>' lowerFcnName '</b>']);
    end
end

% Construct the header.
[pathname filename]=fileparts(topic);
outStr = [makeHeader(filename, '', num2str(count)) sprintf('\n<b>%s</b>',topic)];
% Make "see also", "overloaded methods", etc. hyperlinks.
helpStr = makehelphyper('helpwin', '', fcnName, helpStr);
% Write out the help material
outStr = [outStr sprintf('<br> <br>\n%s\n',helpStr)];


return
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
function h = makeHeader(leftText,centerHtml,rightText)
%--------------------------------------------------------------------------
% Left chunk.
leftData = ['<b>' leftText '</b>'];
leftData = ['<td>&nbsp;' leftData '</td>'];

% Center chunk.
centerData = ['<td valign="left">' centerHtml '</td>'];

% Right chunk.
rightData = ['<b>' rightText '</b>'];
rightData = ['<td align=right>' rightData '&nbsp;</td>'];

beginTable = '<table width="100%" border=0 cellspacing=0 bgcolor=d0d0f0><tr>';
endTable = '</tr></table>';

h = [beginTable leftData centerData rightData endTable];
end
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
function outstr = fixsymbols(instr)
%--------------------------------------------------------------------------
% Convert < and > to html-friendly characters; but do not convert them
% if they are part of an html tag, because then they will show up as
% symbols rather than links inside the Help browser.  So, we need to first
% split the string up so that we don't look at the html tags when replacing
% the symbols with the html-friendly characters.

instr  = strrep(instr, '&', '&amp;');
expr = '<a\shref\s?=\s?.*?>.*?</a>';
[startpos endpos matches] = regexp(instr, expr,'start','end','match');

if ~isempty(startpos)
    startpos = [startpos-1 length(instr)];
    endpos = [1 endpos+1];
    outstr = '';
    for i=1:length(startpos)
        % convert < and > to html-friendly characters
        segment = instr(endpos(i):startpos(i));
        segment = strrep(segment, '<', '&lt;');
        segment = strrep(segment, '>', '&gt;');

        if i<=length(matches)
            % if there are any help links, convert them to helpwin
            outstr  = [outstr segment regexprep(matches{i}, '\<matlab:help\>', 'matlab:helpwin')];
        else
            % append the rest of the help text
            outstr = [outstr segment];
        end
    end
else
    instr  = strrep(instr, '<', '&lt;');
    outstr = strrep(instr, '>', '&gt;');
end
end


% This function copied from MATLAB helptool private folder
%--------------------------------------------------------------------------
function hyperHelp = makehelphyper(actionName, pathname, fcnName, helpStr)
%MAKEHELPHYPER Reformat help output so that the content has hyperlinks

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2007-03-28 15:22:51 $

% Isolate the function name in case a full pathname was passed in
[unused fcnName] = fileparts(fcnName);  %#ok

% Make "see also" references act as hot links
CR = sprintf('\n');
% use sprintf to force translation 
seeAlso = sprintf('See also');
lengthSeeAlso = length(seeAlso);
xrefStart = strfind(helpStr, seeAlso);
% If we are on a Japanese machine but the m-file help is not translated we
% need to look for the English string.
if isempty(xrefStart) && strcmp('See also',seeAlso) == 0
    xrefStart = strfind(helpStr, 'See also');
    % If we found it, reset lengthSeeAlso
    if ~isempty(xrefStart)
        lengthSeeAlso = length('See also');
    end
end
if ~isempty(xrefStart)
    seeAlsoIndex = xrefStart(end);

    % Note : '.' is now valid (for MATLAB class syntax); '>' is valid for subfunctions.
    nameChars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_/.>';  
    delimChars = [ ', ' CR ];
    % Determine start and end of "see also" portion of the help output
    pieceStr = helpStr(seeAlsoIndex+lengthSeeAlso : length(helpStr));
    notePos = min([findstr(pieceStr, sprintf('Overloaded')) findstr(lower(pieceStr), sprintf('note:'))]);
    crPos = regexp(pieceStr,'\n\s*\n');  % blank line with a possible space (or more).
    if isempty(notePos) && isempty(crPos)
        xrefEnd = length(helpStr);
        trailerStr = '';
    elseif ~isempty(notePos)
        xrefEnd = seeAlsoIndex+lengthSeeAlso + notePos(1) - 2;
        trailerStr = pieceStr(notePos(1):length(pieceStr));
    else
        xrefEnd = seeAlsoIndex+lengthSeeAlso + crPos(1) - 1;
        trailerStr = pieceStr(crPos(1):length(pieceStr));
    end

    % Parse the "See Also" portion of help output to isolate function names.
    seealsoStr = '';
    word = '';
    for chx = seeAlsoIndex+lengthSeeAlso : xrefEnd
        if length(findstr(nameChars, helpStr(chx))) == 1
            word = [ word helpStr(chx)];
        elseif (length(findstr(delimChars, helpStr(chx))) == 1)
            if length(word) > 0 
                % This word appears to be a function name.
                % Make link in corresponding "see also" string.
                fname = word;
                if strcmp(upper(word), word) == 1
                    % only lowercase this if it's all uppercase.
                    fname = lower(word);
                end
                suffix = '';
                if fname(length(fname)) == '.'
                    % Don't hyperlink the last period of the word.
                    fname = fname(1:length(fname)-1);
                    suffix = '.';
                end
                % Make sure the function exists before hyperlinking it.
                if strcmp(word,'and') == 0 
                    % don't even try to hyperlink 'and'... otherwise, we
                    % will hyperlink all words in See Also.
                    seealsoStr = [seealsoStr '<a href="matlab:' actionName ' ' fname '">' fname '</a>' suffix];
                else
                    seealsoStr = [seealsoStr word];
                end
            end

            
            seealsoStr = [seealsoStr helpStr(chx)];
            word = '';
        else
            seealsoStr = [seealsoStr word helpStr(chx)];
            word = '';
        end
    end
    % Replace "See Also" section with modified string (with links)
    helpStr = [helpStr(1:seeAlsoIndex+lengthSeeAlso -1) seealsoStr trailerStr];
end

% If there is a list of overloaded methods, make these act as links.
overloadPos =  findstr(helpStr, xlate('Overloaded'));
if strcmp(actionName,'doc') == 1
    textToFind = ' doc ';
    len = 5;
else
    textToFind = ' help ';
    len = 6;
end

if length(overloadPos) > 0
    pieceStr = helpStr(overloadPos(1) : length(helpStr));
    % Parse the "Overload methods" section to isolate strings of the form "help DIRNAME/METHOD"
    overloadStr = '';
    linebrkPos = find(pieceStr == CR);
    lineStrt = 1;
    for lx = 1 : length(linebrkPos)
        lineEnd = linebrkPos(lx);
        curLine = pieceStr(lineStrt : lineEnd);
        methodStartPos = findstr(curLine, textToFind);
        methodEndPos = length(curLine) - 2;
        if (length(methodStartPos) > 0 ) && (length(methodEndPos) > 0 )
            linkTag = ['<a href="matlab:' actionName ' ' curLine(methodStartPos(1)+len:methodEndPos(1)+1) '">'];
            overloadStr = [overloadStr curLine(1:methodStartPos(1)) linkTag curLine(methodStartPos(1)+1:methodEndPos(1)+1) '</a>' curLine(methodEndPos(1)+2:length(curLine))];
        else
            overloadStr = [overloadStr curLine];
        end
        lineStrt = lineEnd + 1;
    end
    % Replace "Overloaded methods" section with modified string (with links)
    helpStr = [helpStr(1:overloadPos(1)-1) overloadStr];
end

% If this topic is a Contents.m file, scan it for function lists, and
% modify function names to act as active links.
if (strcmpi(fcnName, 'Contents')) || (length(findstr(helpStr,'is both a directory and a function'))) || (strcmp(fcnName,'simulink')) || (strcmp(fcnName,'debug'))
    TAB = sprintf('\t');
    helpStr = strrep(helpStr, TAB, ' ');
    modHelpStr = '';
    linebrkPos = find(helpStr == CR);
    lineStrt = 1;
    for lx = 1 : length(linebrkPos)
        lineEnd = linebrkPos(lx);
        curLine = helpStr(lineStrt : lineEnd);
        hyphPos = findstr(curLine, ' - ');
        if any(hyphPos)
            nonblankPos = find(curLine ~= ' ');
            if curLine(nonblankPos(1)) == '-'
                modHelpStr = [modHelpStr curLine];
            else
                % start with the hyphen position, and walk backwards to the
                % first nonblank position.
                for i = hyphPos(1):-1:nonblankPos(1)
                    if (curLine(i) ~= ' ') && (curLine(i) ~= ','), break, end;
                end
                fname = curLine(nonblankPos(1):i);
                fname = strrep(fname, '&amp;', '&');
                remainder = curLine(i+1:end);
                try
                    % If there is any help for this name, insert a link for it.
                    % First determine if we need to qualify the name to put
                    % in the link by testing whether or not the function
                    % exists under the same directory.
                    qualified_name = fname;
                    fullpath = which(fname);
                    if ~isempty(pathname) && ~isempty(fullpath)
                        startpos = strfind(fullpath, [pathname filesep]);
                        if isempty(startpos)
                            % Maybe the function is shadowed... if so, we
                            % need to add a prefix to the help command so
                            % we get the correct help.
                            fullpath = which('-all',fname);
                            if length(fullpath)>1
                                for entry=2:length(fullpath)
                                    pathentry = fullpath{entry};
                                    startpos = strfind(pathentry, [filesep pathname filesep]);
                                    if ~isempty(startpos)
                                        pname = strrep(pathentry(startpos+1:end),'\','/');
                                        [pname unused] = fileparts(pname);  %#ok
                                        if ~isempty(pname)
                                            qualified_name = [pname '/' fname];
                                        end
                                    end
                                end
                            end
                        end
                    end
                        
                    % If this is a valid function name, convert it to a
                    % hyperlink.
                    if isHyperlinkable(fname)
                        fnameLink = strrep(qualified_name,'''','''''');
                        if remainder(1) == ','
                            % Sometimes there are two names separated by a comma.
                            [fname2, remainder2] = strtok(remainder(3:end));
                            fnameLink2 = strrep([pathname '/' fname2],'''','''''');
                            modHelpStr = [modHelpStr curLine(1:nonblankPos(1)-1) '<a href="' 'matlab:' actionName ' ''' fnameLink '''">' fname '</a>, ' '<a href="' 'matlab:' actionName ' ''' fnameLink2 '''">' fname2 '</a>' remainder2];
                        else
                            modHelpStr = [modHelpStr curLine(1:nonblankPos(1)-1) '<a href="' 'matlab:' actionName ' ''' fnameLink '''">' fname '</a>' remainder];
                        end
                    else
                        modHelpStr = [modHelpStr curLine];
                    end
                catch
                    % Just in case an error occurred during the helpfunc call, don't try to
                    % hyperlink anything.
                    modHelpStr = [modHelpStr curLine];
                end
            end
        else
            modHelpStr = [modHelpStr curLine];
        end
        lineStrt = lineEnd + 1;
    end
    helpStr = modHelpStr;
end
hyperHelp = helpStr;
end

function shouldLink = isHyperlinkable(fname)
% Make sure the function exists before hyperlinking it.
shouldLink = 0;
fnameType = exist(fname);  %#ok
if fnameType == 2 || fnameType == 3 || fnameType == 5 || fnameType == 6 || fnameType == 7
    % help may exist for M-files, MEX-files, P-files (which should have a corresponding
    % M-file with help), built-ins, and directories.  But not variables, MDL-files, or Java classes.
    shouldLink = 1;
elseif fnameType == 8
    % this means the file is either a UDD class or a Java class
    fullpath = which(fname);
    fType = exist(fullpath,'file');
    if fType == 2 || fType == 3 || fType == 7
        % Since it's an M-file, P-file, or MEX-file, it should have some help.
        shouldLink = 1;
    end
elseif fnameType == 0
    % Check to see if it's a UDD method
    dotpos = findstr(fname,'.');
    slashpos = findstr(fname,'/');
    if length(dotpos) == 2 || (length(dotpos) == 1 && length(slashpos) == 1 && dotpos(1)<slashpos(1))
        fullname = strrep(fname, '.', '/');
        fullname = ['@' fullname(1:dotpos(1)) '@' fullname(dotpos(1)+1:end)];
        if ~isempty(which(fullname))
            shouldLink = 1;
        end
    end
end
end





