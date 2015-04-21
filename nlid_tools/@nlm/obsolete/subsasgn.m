function Aout=subsasgn(A,S,B)
% subsasgn for nlm objects
% B must be a system or cell array of systems;

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

[nx(1),nx(2)]=size(B);
for i=1:2,
  ix{i}=1:nx(i);
end

for i=1:length(S.subs),
  if ischar(S.subs{i}),
    ix{i}=1:nx(i);
  else
    ix{i}=S.subs{i};
  end
end
Aout=A;
E=A.Elements;   
if strcmp (S.type, '()')
  % RHS must be cell array
  E (ix{1},ix{2})=B;
  Aout.Elements= { E };
  
elseif strcmp(S.type,'{}')
  E (ix{1},ix{2})={B };
  Aout.Elements=  E ;
else
  error ('format not defined for nlm objects');
  
end
% @nlm/subsasgn
