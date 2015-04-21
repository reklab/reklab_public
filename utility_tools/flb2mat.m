function D= flb2mat_m (fname, option, caseNum)
% FLBTOMAT_M.M ...
%
%  read data from a  FILELB file
%%   Usage:  D = flbtomat(file, option, case)
%  outputs is a strucutre containing details and daa
%
% % Usage
%  inputs are:
%          file - name of filelb file
%          option  =  'index' - generate index of flbfile
%                     'length' - number of cases in file
%                     'read_case' - read one or more cases
%                     'read_all' - read all cases
%          caseNum = caseNumber to read
%          option  = q for quit mode index
%
%  Calling flbtomat with no parameters gives a help message
%
%
% Example
%
% Read case umber 3 from flb file test.flb
%
%  D = flb2mat ('test.flb','read_case', 3);
%
% Reall all data into a strucutre array
%
%  D = flb2mat('test.flb','read_all')
%
% Get idex of flbfile
%
% idex = flb2mat('test.flb','index');
%
% Determine number of cases
%
% nCase = flb2mat('test.flb','length');

%
% See also
%
%% AUTHOR: Robert Kearney
%% $DATE: 01-Feb-2006 17:36:09 $
%% $Revision: 1.3 $
%% DEVELOPED: 7.1.0.246 (R14) Service Pack 3
%% FILENAME: flbtomat_m.m.m
%Copyright 2002-2005 McGill University This file is part of the CMB Proteomics PipeLine

%% open file
fid=fopen(fname);
if fid < 1,
    error('error opening file');
end
if nargin==1,
    option='index';
end

%%
switch lower(option)
    case 'index'
        D=flbio(fid,'read_index');
        if nargout==0,
            struct_disp(D);
        end
    case 'length'
        tempD=flbio(fid,'read_index');
        D=length(tempD);
    case 'read_case'
        for iCase=1:caseNum-1,
            tempD=flbio(fid,'skip_case');
        end
        D=flbio(fid,'read_case');
        
    case 'read_all'
        D=flbio(fid,'read_all');
    otherwise
        error([' mat2flb option not defined:' option]);
end
fclose(fid);



return
