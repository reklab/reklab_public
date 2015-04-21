function [Rval] = pnv(sys,mode, case_flag)
% pnv - returns public properties or vlaues of object sys. 
% sys - object
% mode - mode to return [names/values]
% case - 'lower' returns values as lower case 

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if nargin < 3,
    case_flag='lower';
end
%
% Customize for particular application here
%
parent=sys.nltop; % Parent object
Npublic = 8;  % Number of MS specific public properties
% Define admissible values here 
    AsgnVals = {'cell array of strings'; ...
      'cell array of strings'; ...
      'real [1]'; ...
      'real[x]'; ...
      'real [0]'; ...
      'numeric array'; ...
      'numeric array[x]' ; ...
      'size of data set' };




switch lower(mode)
case 'names'
    
    
    %
    %   [PROPS,ASGNVALS] = pnv (SYS,'names', 'true')  returns the list PROPS of
    %   public properties of the object SYS (a cell vector), PROPS contains the true case-sensitive property names.
    %   These include the public properties of SYS's parent(s).
    %
    %   [PROPS,ASGNVALS] = PNAMES(SYS,'lower')  returns lowercase property
    %   names.  This helps speed up name matching in GET and SET.
    %
    % Define propery names here
    
    Props = fieldnames(sys);
    Props = Props(1:length(Props)-1); 
    
    
    ;
    % Get parent properties and values
    
    Rval = [Props ; pnv(parent,'names')];
    if strcmp(lower(case_flag),'lower'),
        Rval = lower (Rval);
    end
case 'avalues'
    
    %   Return the assignable values ASGNVALS for public properties (a cell vector 
    %   of strings) of hte object and its parent
    
    
    % Add parent  admissible values
    ParentAvals = pnv(parent,'avalues');
    Rval = [AsgnVals ; ParentAvals];
    
case 'values'
    %   VALUES = pnv(SYS,'values')  returns the list of values of all
    %   public properties of the object SYS.  VALUES is a cell vector.
    %
    
    
    
    % Values of public nldat specific properties
    Values = struct2cell(sys);
    Values = Values(1:Npublic);
    
    % Add parent properties
    Rval = [Values ; pnv(parent,'values')];
otherwise
    error ('option not defined for pnv');
end



% end cor/pnv.m
