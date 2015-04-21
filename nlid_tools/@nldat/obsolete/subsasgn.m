function Aout=subsasgn(A,S,B)
% subsasgn for nldat objects
% V01-01 12 Nov 98
%

% Copyright 1998-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

B=nldat(B);
if strcmp (S.type, '()')
  [nx(1),nx(2),nx(3)]=size(B.Data);
  for i=1:3,
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
  Aout.Data(ix{1},ix{2},ix{3})=B.Data;
  %
  % Make sure output channel names are defined
  cnames = B.ChanNames;;
  l1=min(ix{2});
  l2=max(ix{2});
  for l=l1:l2,
    cn1=cat(2,'x',int2str(l));
    cnames=cat(2,cnames,{cn1}); 
  end
  Aout.ChanNames=cnames;
  
  
elseif strcmp(S.type,'{}')
  error ('{} format not defined for nldat objects');
  
end
