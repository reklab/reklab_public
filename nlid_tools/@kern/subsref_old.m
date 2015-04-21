function b = subsref (A,S);
% syntax  b = subsref (A,S);

% Copyright 1991-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 


if strcmp (S.type, '()'),
    K=get(A,'data');
   [nx(1),nx(2),nx(3)]=size(K);
   for i=1:3,
      ix{i}=1:nx(i);
   end
   for i = 1:length(S.subs),
      if ischar(S.subs{i}),
         ix{i}=1:nx(i);
      else
         ix{i}=S.subs{i};
      end
    end
    b=A;
    set(b,'data',K(ix{1},ix{2},ix{3}));
 
 elseif strcmp(S.type, '{}'),
   ix=S.subs{1};
   b=double(K{ix+1});
elseif strcmp(S.type,'.')
    b=A.(S.subs);
else
   b=NaN;
   error ('subsref function not yet implemented');
end
return
