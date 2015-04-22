function error_flag = arg_parse_c (  caseMode, options, varargin)
% arg_parse ( caseMode, options, varargin) -Parse and decode variable number of inputs with case control
% o
% Usage
% caseMode - detemrines have variable names are handle
%   'exact' - varargin must match options exactly; case preserved;
%   'lower' - case ignored in matching varargin to options; variable names returne in lower case
%   'upper - case ignored in matching varargin to options;variable names returne in upper case
% options cell array of form:
%       { {varname1 default1} {option2 default2} ... }
% varargin - name value pairs
%
% arg_parse - matches names in options to thise in varargin.
%               if option is in varargin then the varible is set to
%               specified value;
%               if option is not in varargin then varialbe is set to
%               default value specified in options
%               if varargin name is not in options an error is generated.
%
% Example
%
% for options = { 'var1' 1} { 'var2' 'test'} }
%
%   arg_parse_c ( 'lower' , options )
%       sets var1=1 and var2='test' in calling routine
%   arg_parse_c ( 'lower' , options, 'var1', 10,'var2','new string'  )
%        sets var1=1 and var2='test' in calling routine
%
% See also
%
%% AUTHOR: Robert Kearney
%% $DATE: 23-Aug-2006 11:24:30 $
%% $Revision: 1.5 $
%% DEVELOPED: 7.2.0.232 (R2006a)
%% FILENAME: arg_parse_c.m

% Set default values
error_flag=0;
nopt=length(options);
for i=1:nopt,
    names{i} = getName(caseMode,(char(options{i}(1,1))));
    defaults(i)=options{i}(1,2);
    assignin ('caller',char(names{i}), defaults{i});
end
if (nargin <3),
    return
elseif isempty(varargin)
    return
end

% Take care of nested calls with variable number of parameters
% where args appear as one element of a cell array

while length(varargin)==1 & iscell (varargin{1}),
    varargin=varargin{1};
end
narg=length(varargin);
if length(varargin) > 0 & strcmp('?',varargin(1)),
    s=dbstack;
    arg_help(s(1).name, options);
    error_flag=1;
    return
end

%
% Parse variable arguments
%
for i=1:2:narg,
    choice=getName(caseMode,(char(varargin{i})));
    j=find(strncmp(choice, names, length(choice)));
    if j > 0,
        assignin('caller',choice,varargin{i+1});
    else
        disp (['Invalid option:' choice]);
        vararg_help(options);
        error(['Bad Option:' choice]);
    end
end
end

function vararg_help (options)
fprintf (1, '\nOptions are:')
for i=1:length(options),
    fprintf(1,'\n\t');
    fprintf(1,'%c',char(options{i}(1,1)))
end
fprintf(1,'\n');
return
end


%%% function getName - get name with appropriae case conversion
function nOut =getName (caseMode,nIn);
switch lower(caseMode),
    case 'lower'
        nOut=lower(nLin);
    case 'upper'
        nOut=upper(nIn);
    otherwise
        nOut=nIn;
end
end









% Created with NEWFCN.m by Frank González-Morphy
% Contact...: frank.gonzalez-morphy@mathworks.de
% ===== EOF ====== [arg_parse_c.m] ======
