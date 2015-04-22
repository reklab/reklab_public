function p = nlident (p, z,  varargin)
% CONSTRUCT an pdf  function object
%
% Setup default values
%
% $Revision: 1.9 $
% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

if nargin < 2,
  disp('NLIDtakes two inputs for pdf objects: pdf, z' );
elseif nargin > 2,
  set(p,varargin);
end
assign (p.Parameters); 
tempComment=p.Comment; 
if isa(z,'nldat') | isa (z,'double'),
  if isa(z,'double'),
    z=nldat(z);
  end
  
  x=double(z);


  % Number of bins

  if isnan (NBins)
    NBins=length(x)/10.;
    set(p,'NBins',NBins);
  end
  if isnan(BinMin),
    BinMin=double(min(z.Data));
    set(p,'BinMin',BinMin);
  end
  if isnan(BinMax),
    BinMax=double(max(z.Data));
    set(p,'BinMax',BinMax);
  end
  if BinMin >=BinMax,
      BinMax=BinMin+1;
     
  end

  % 
  pdfType=lower(pdfType); 
  [dist,dom]=pdfnl(x,NBins,pdfType,BinMin,BinMax);
  incr=dom(2)-dom(1);
  switch pdfType
      case 'cumulativeprobability'
          set(p,'Data',dist','ChanNames',{'CumulativeProbabilty'}, ...
	  'ChanUnits','Probabilty','Comment','Cumultive Probabilty');
      case 'density'
      set(p,'Data',dist','ChanNames',{'Density'}, ...
	  'ChanUnits','count','Comment','Density');
    case 'frequency'
      set(p,'Data',dist','ChanNames',{'Frequency'}, ...
	  'ChanUnits','count','Comment','Frequency');
    case 'probability'
      set(p,'Data',dist','ChanNames',{'Probabilty'}, ...
	  'ChanUnits',' ','Comment','Probability');
   
  end
  set (p,'DomainName',[ 'Value of ' char(z.ChanNames)]);
  set (p,'DomainStart',min(dom), 'DomainIncr',incr);
  set(p,'DomainValues',dom);

  %set(p,'NBins',NBins,'BinMin',BinMin,'BinMax',BinMax);
  %
  %
elseif isa(z,'randv'),
  p=randv2pdf(z,pdfType, NBins, BinMin,BinMax);
  
end
if ~strcmp(tempComment,p.Comment),
    set(p,'Comment',tempComment);
end
