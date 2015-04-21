function status = mat2flb(fileName, option, D, varargin )
% MAT2FLB ...
% function status = mat2flb (fileName, option, D )
% Store matlab files as  filelb files
% Usage
%

%          fileNamne - name of filelb file
%          option - 'append' to current file
%                 - 'overwrite' data in current file
%                 - 'newreal' - add a realization to last case
%          D-  data to store with a matrix of form  D(nSamp, nChan, NReal)
%               or a data strucutre with format
%                D.comment - string variable
%                D.chanName - cell array of channel names
%                D.domainIncr - double
%                D.domainName - string
%                D.domainStart - double
%                D.Data - matrix of format (nSamp, nChan, nReal);
%
%          varargin - name value pairs withthe following possibilites
%          comment string
%          chanNames  cell array of strings
%          domainIncr - double
%          domainName - string
%          domainStart - double
%
%  Only  the filename and matrix paramters are mandatory
%  HELP is dislayed if flbtomat is called with no parameters
%   ...

%
% Example
%
% Store the strucutre D in file test.flb overwriting previous contents
%
%   mat2flb ('test.flb', 'overwrite' , D);
%
% Append the strucutre D to contensts of test.flb
%
%   mat2flb('test.flb', 'append', D);
%
% Add a realization to the current case
%
%   mat2flb('test.flb', 'append', D);
%   mat2flb('test.flb', 'newreal', D);
%
%
% See also
%
%% AUTHOR: Robert Kearney
%% $DATE: 02-Feb-2006 10:46:49 $
%% $Revision: 1.3 $
%% DEVELOPED: 7.1.0.246 (R14) Service Pack 3
%% FILENAME: mat2flb.m
%Copyright 2002-2005 McGill University This file is part of the CMB Proteomics PipeLine
%% set default channel names

% DL made changes to line 108
if ~isstruct(D),
    dataSize=size(D);
    nSamp=dataSize(1);
    nChan=dataSize(2);
    if nChan > nSamp &  nChan > 10,
        warning (' nCHan > nSamp. Transposing');
        D=D';
        nChan=nSamp;
    end
else
    dataSize=size(D.Data);
    nChan=dataSize(2);
end

%% Set default values

for iChan=1:nChan,
    chanName{iChan}= [ 'x' int2str(iChan)];
end


argOptions = { { 'comment' 'default comment' } ...
    { 'domainname' 'default domain' } ...
    { 'domainstart' 0 } ...
    { 'domainincr' 1 } ...
    { 'channame' chanName}};

arg_parse(argOptions, varargin);


%% Check input format and set default values
% Input is a matrix so generate data structure using either default values
%% or those specified in input
if ~isstruct(D),
    xtemp=D;
    clear D;
    % DL changed from D.Data=Dtemp;
    D.Data=xtemp;
end
%% Input is a structure so check structure format and set missing elements
%% to default values
if isstruct(D),
    % check for data
    if ~isfield (D,'Data');
        error('Data field is missing. Fatal error');
    end
    % comment field
    if ~isfield (D, 'comment'),
        D.comment= comment;
    else
        dataSize=size(D.Data);
        nChan=dataSize(2);
    end

    if ~isfield (D,'domainIncr'),
        D.domainIncr= domainincr;
    end
    if ~isfield (D,'domainName'),
        D.domainName= domainname;
    end
    if ~isfield(D,'domainStart'),
        D.domainStart=domainstart;
    end
    if ~isfield(D,'chanName'),
        D.chanName=channame
    elseif length(D.chanName) ~=nChan,
        error('Length chanName field does not match number of channels in data');
    end
end

switch lower(option)
    case 'append'
        fid=fopen(fileName, 'a');
        if fid<1,
            error('Error opening file');
        end
        s=flbio(fid,'write_case', D);
    case 'overwrite'
        disp('This will delete the current file');
        if (input_l('Procede', false)),
            fid= fopen(fileName,'w+');
            if fid<1,
                error('Error opening file');
            end
        else
            disp('operation terminated');
            return
        end
        s=flbio(fid,'write_case', D);
    case 'newreal'
        fid= fopen(fileName,'r+');
        if fid<1,
            error('Error opening file');
        end
       s =flbio(fid,'new_real', D); % Add a realization 
    otherwise
        error([' mat2flb option not defined:' option]);
end
fclose(fid);






