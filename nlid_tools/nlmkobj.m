function O = nlmkobj ( obj, varargin )
% NLMKOBJ - Standard call for NLID objects
% Parse input
%
% $Revision: 1.5 $

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt

objclass=class (obj);
O=obj;
nArg=nargin;
% Only one argument so simply return it
if nargin == 1 | isempty(varargin),
    return
end
z=varargin{1};

if ischar(z),
    % Second Argument is character so set  values  
    set(O,varargin{:});
    
else 
    % Do identification
    if nargin==2,
        O=nlident (O,z);
    else
        O=nlident (O,z, varargin{2:end});
    end
    % z is data so do identification
end
end

% nlmkobj

