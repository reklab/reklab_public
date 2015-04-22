function [s,m]= rsync(passWord, optionList, source, dest)
% RSYNC - matlab wrapper for rsyn  
%  
%   ... 
% Usage 
% Input paramters:
%   password - rsync password
%   option list - rsync options
%   source - source directory
%   dest - destinaton directory
% Returns:
%   s - stats 0=success otherwise error status
%   m = message

% 
% Example 
% 
% See also 
% 
%% AUTHOR: Robert Kearney 
%% $DATE: 20-Oct-2006 13:56:05 $ 
%% $Revision: 1.1 $ 
%% DEVELOPED: 7.3.0.267 (R2006b) 
%% FILENAME: rsync.m 
%Copyright 2002-2005 McGill University This file is part of the CMB
%Proteomics PipeLine

switch computer
    case 'PCWIN'
        rsyncCMD = 'rsync';
        % change windows d:/ to /cygdrive/d/
          source=strrep(source,'\','/'); 
          if strcmp(source(2),':'),
            source=strrep(source,':','');          
            source = [ '/cygdrive/' source];
            
        end
    otherwise
        error('system not yet supported');
end
space=' ';

cmd= [rsyncCMD   space  optionList  space quote(source)  space dest];

%%
setenv('RSYNC_PASSWORD',passWord);
setenv('CWRSYNCHOME','C:\PROGRAM FILES\CWRSYNC');
setenv('CYGWIN','nontsec')
setenv('HOME', [getenv('HOMEDRIVE') getenv('HOMEPATH')]);
oldPath=getenv('PATH');
setenv('PATH', [ 'C:\PROGRAM FILES\CWRSYNC\BIN;' getenv('PATH')]);

[s,m]=system(cmd);

setenv('PATH',oldPath);


    




% Created with NEWFCN.m by Frank González-Morphy  
% Contact...: frank.gonzalez-morphy@mathworks.de  
% ===== EOF ====== [rsync.m] ======  
