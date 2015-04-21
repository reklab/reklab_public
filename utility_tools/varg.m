function [N,V] = varg (N, V )
% examine varargin and correct for nesting calls 
%

if N >1
   while isa(V{1},'cell'),
      V=V{1};
   end
   N=length(V)+1;
end

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 
