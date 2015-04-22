function s1 = deblankl(s)
%DEBLANKL Remove leading blanks.
%   DEBLANKL(S) removes leading blanks from string S.
%
%   DEBLANKL(C), when C is a cell array of strings, removes the leading
%   blanks from each element of C.


%   $Revision: 1.1 $  $Date: 2004-01-20 02:38:40 $

% The cell array implementation is in @cell/deblank.m

if isempty(s)
   s1 = s([]);
else
   if ~isstr(s),
      warning('MATLAB:deblankl:NonStringInput','Input must be a string.')
   end
   
   % Remove trailing blanks
   [r,c] = find( (s~=0) & ~isspace(s) );
   if isempty(c),
      s1 = s([]);
   else
      s1 = s(:,1:min(c));
   end
end
