function [index,status] = pmatch(StrCell,str)
%PMATCH  Case insensitive string matching with automatic completion.
%
%   [INDICES,STATUS] = PMATCH(STRCELL,STR)  takes a cell array 
%   STRCELL containing reference names and a string STR to be 
%   matched.  It returns an index vector INDEX such that 
%   STRCELL(INDEX) are all the strings matching STR.
%   
%   STATUS is the empty matrix if exactly one cell matches each
%   string.  Otherwise it returns a "standard" error message.
%
%   Note:  PMATCH is optimized for the LTI object and uses only
%   the first two characters of STR for comparison purposes.

%   Author(s): A. Potvin, 3-1-94
%   Revised: P. Gahinet, 4-1-96
%   Copyright (c) 1986-98 by The MathWorks, Inc.
%   $Revision: 1.2 $

%   Note: replace nchars = min(l,2) by nchars = l below to obtain 
%   following general behavior:
%        StrCell = {'CurrentMenu'; 'Color'; 'CurrentObject'}
%        STR = 'current' 
%           INDICES = []
%           STATUS  = 'Ambiguous property name: supply more characters'
%        STR = 'currentmenu' 
%           INDICES = 1
%           STATUS  = []
%        STR = 'foobar' 
%           INDICES = []
%           STATUS  = 'Invalid object property'


index = []; 
status = [];
if ~isstr(str),
   status = 'Property names must be strings';
   return
elseif ndims(str)>2 | min(size(str))>1,
   status = 'Property names must be single-line strings';
   return
end

% STR is a single-line string
str = str(:)';
s = lower(str);
ls = length(s);

% Set number of chars used to identify name
nchars = min(ls,10);  

% Find all matches in StrCell
% RE: Assumes all strings in StrCell are lower case
iTemp=strcmp(StrCell,str);
index=find(iTemp); 


% Error handling 
nhits = length(index);
if nhits==0,
   status = ['Invalid object property "' str '"'];
   index = [];  return
elseif nhits>1,
   status = ['Ambiguous property name "' str '". Supply more characters.'];
   index = [];  return
end


% end pmatch.m
