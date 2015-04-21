function selectString = guiInput_List(promptString, listString)
% GUIINPUT_LIST - select from a strong from a list
%  
%   ... 
% Usage 
%
% selectString = guiInput_List(promptString, listString)
% 
% Example 
% 
% See also 
% 
%% AUTHOR: Robert Kearney 
%% $DATE: 22-Nov-2007 10:43:46 $ 
%% $Revision: 1.2 $ 
%% DEVELOPED: 7.5.0.342 (R2007b) 
%% FILENAME: guiInput_List.m 
%Copyright 2002-2005 McGill University This file is part of the CMB Proteomics PipeLine 

listLen=length(listString);
listSize=[160 listLen*20];
drawnow
[select, OK] = listdlg ( 'PromptString', promptString, 'ListString' ,listString, ...
    'SelectionMode', 'single', 'Name', 'guiInput_List', 'ListSize', listSize);
drawnow
if OK,
    selectString= listString{select};
else
    error ('Cancel');
end






% Created with NEWFCN.m by Frank González-Morphy  
% Contact...: frank.gonzalez-morphy@mathworks.de  
% ===== EOF ====== [guiInput_List.m] ======  
