function e = subsref (N,S);
% overlaid subsref for nlm objects
% N(i,j) returns system description ith object in jth path

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

E=N.Elements;
[nx(1),nx(2)]=size(E);
for i=1:2,
  ix{i}=1:nx(i);
end
for i = 1:length(S.subs),
  if ischar(S.subs{i}),
    ix{i}=1:nx(i);
  else
    ix{i}=S.subs{i};
  end
end
if strcmp (S.type, '()'),
  % Return Cell
  e=E(ix{1},ix{2});
elseif strcmp (S.type, '{}'),
  % Return Contents
  e=E{ix{1},ix{2}};
else
  e=NaN;
  error ('NLM subsref not yet implemented for this syntax');
end
return
% nlm/subsref
