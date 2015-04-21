function datout = subsref (datin,S);
% V01-02 REK 22 Oct 99 Fix problem when domain values are specifed
% V01-03 REK 20 Dec 02 Add support for '.' format
% $Revision: 1.3 $
% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

while length(S) >1,
    datin=subsref(datin,S(1));
    Snew=S(2:length(S));
    S=Snew;
end
datout=datin;
if strcmp (S.type, '()'),
   [nx(1),nx(2),nx(3)]=size(datin.Data);
   for i=1:3,
      ix{i}=[1:nx(i)];
   end
   for i = 1:length(S.subs),
      if ischar(S.subs{i}),
         ix{i}=1:nx(i);
      else
         ix{i}=S.subs{i};
      end
    end
    b=datin;
    d=datin.Data; 
    dout=d(ix{1},ix{2},ix{3});
    set(datout,'Data',dout);
    % Upddate Domain Values
    if ~isnan (datin.DomainValues),
       adv=datin.DomainValues;
       datout.DomainValues=adv(ix{1});
       datout.DomainStart=datout.DomainValues(1);
    else
       datout.DomainStart=datin.DomainStart + (ix{1}(1)-1)*datin.DomainIncr;
    end
  
    acn=datin.ChanNames;
    bcn=acn(ix{2});
    datout.ChanNames=bcn;
elseif strcmp(S.type,'.')
    datout=datin.(S.subs);
   
else
   datout=NaN;
   error ('subsref function not yet implemented');
end
return
