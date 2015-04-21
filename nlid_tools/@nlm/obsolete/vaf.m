function vf = vaf (model,data,spare)
% VAF function for NLM class models
%
% vf = vaf(model,data)
%
% where model is any model of class nlm
% and data is an nldat object containing [input,output]

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if nargin == 1
  error('I/O data must be supplied if first argument is of class NLM');
end

if nargin == 2
  % data contains both input and output
  if isa(data,'nldat');
    ddata = double(data);
    Ts = data.DomainIncr;
  elseif isa(data,'double')
    ddata = data;
    Ts = 0;
  else
    error('second argument must be either double or nldat');
  end
  [N,m] = size(ddata);
  if min(N,m)==1
    error('must supply both input and output to compute vaf');
  elseif m > 2
    error('To test mimo systems, vaf(model,inputs,outputs)');
  end
  u = nldat(ddata(:,1));
  y = nldat(ddata(:,2));
  if Ts > 0
    set(u,'DomainIncr',Ts);
    set(y,'DomainIncr',Ts);
  end
end

if nargin == 3
  if isa(data,'nldat');
    u = data;
  elseif isa(data,'double');
    u = nldat(data);
    Ts = 0;
  else 
    error('second argument must be either double or nldat');
  end
  if isa(spare,'nldat');
    y = spare;
    if Ts == 0
      Ts = get(y,'domainincr');
      set(u,'domainincr',Ts);
    end
  elseif isa(data,'double');
    y = nldat(spare);
    if Ts > 0
      set(y,'domainincr',Ts);
    end
  else 
    error('third argument must be either double or nldat');
  end
end

yest = nlsim(model,u);
if Ts == 0
  Ts = get(yest,'domainincr');
  set(y,'domainincr',Ts);
end

    
vf = vaf(y,yest);
