function Rval=pnvmain( sys,mode, case_flag, parent, Npublic, AsgnVals,top_flag)
%   [PROPS,ASGNVALS] = pnv (SYS,'names', 'true')  returns the list PROPS of
%   public properties of the object SYS (a cell vector), 
%   PROPS contains the true case-sensitive property names.
%   These include the public properties of SYS's parent(s).
%
%   [PROPS,ASGNVALS] = PNAMES(SYS,'lower')  returns lowercase property
%   names.  This helps speed up name matching in GET and SET.
%


% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


 if nargin <7,
        top_flag=0;
    end
    
switch lower(mode)
   
case 'names'
    
    %
    % Define propery names here
    
    Props = fieldnames(sys);
    if ~top_flag,
        Props = Props(1:length(Props)-1); 
    end
    
    
    % Get parent properties and values
    
    if top_flag,
        Rval= Props;
    else
        Rval = [Props ; pnv(parent,'names')];
    end
    if strcmp(lower(case_flag),'lower'),
        Rval = lower (Rval);
    end
case 'avalues'
    
    %   Return the assignable values ASGNVALS for public properties (a cell vector 
    %   of strings) of hte object and its parent
    % Add parent  admissible values
    if top_flag
        Rval=(AsgnVals);
    else
        
        ParentAvals = pnv(parent,'avalues');
        Rval = [AsgnVals ; ParentAvals];
    end
    
    
case 'values'
    %   VALUES = pnv(SYS,'values')  returns the list of values of all
    %   public properties of the object SYS.  VALUES is a cell vector.
    %      
    % Values of public nldat specific properties
    if (Npublic >0),
        Values = struct2cell(sys);
        Values = Values(1:Npublic);
    else
        Values={};
    end
    
    % Add parent properties
    if top_flag
        Rval=Values;
    else
        Rval = [Values ; pnv(parent,'values')];
    end
    
otherwise
    error ('option not defined for pnv');
end
return



% end fresp/pnv.m
