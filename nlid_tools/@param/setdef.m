function Pout = setdefault (Pin, varargin );
% set the  default value for elements within a parameter set
% Pin - input paramter array
% varagin - name/value pairs to set

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

Pout=Pin;
% Take care of nested calls with variable number of parameters
% where args appear as one element of a cell array
   
while length(varargin)==1 & iscell (varargin{1}),
  varargin=varargin{1};
end

for i=1:2:length(varargin),
   j=pindex(Pin,varargin{i});
   Pout.Default{j}=varargin{i+1};
end



function i = pindex (P,name)
for i=1:length(P),
   pname=P.Name{i};
   if strcmp (pname,name),
      return
   end
end
error (['index' name 'not found']);
