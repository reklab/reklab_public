function vs = nlm2vs(nlmodel, vsin)
% generate Volterra kernels for a NLM model

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 



subsys = get(nlmodel,'elements');
[nparallel,nseries] = size(subsys);
Ts = 0;
thisTs = 0;


% 1.  Figure out the maximum kernel order and memory length
nLags = 0;
Q = 0;
for i = 1:nparallel
  memory = 0;
  order = 1;
  for j = 1:nseries
    ss = subsys{i,j};
    if isa(ss,'irf')
      memory = memory + get(ss,'nlags');
      thisTs = get(ss,'domainIncr');
    elseif isa(ss,'kern')
      memory = memory + get(ss,'nlags');
      order = order * get(ss,'order');
      thisTs = get(ss,'domainIncr');
    elseif isa(ss,'polynom')
      order = order * get(ss,'order');
    else 
      error('nlm object contains unconvertable element');
    end
    if Ts == 0
      Ts = thisTs;
    elseif Ts~=thisTs
      error('Elements have different sampling rates');
    end
  end
  nLags = max(nLags,memory);
  Q = max(Q,order);
end

% compute number of kernel elements (since we store the whole kernels, we
% can't make use of symmetries).
nLags = max(1,nLags);
M = 1;
for i = 1:Q
  M = M + nLags^i;
end

if ((Q > 3)|(M>10000))
  % warn the user before proceeding as things could get big and ugly!
  disp(['This will create kernels of order 0 through ' num2str(Q)]);
  disp(['with memory lengths of ' num2str(nLags) ' samples']);
  proceed = input('Continue ([y]/n)','s');
  if isempty(proceed)
    proceed = 'y';
  end
  if proceed(1) == 'n'
    return
  end
end

% 1a allocate storage for the kernels.
kernels = cell(Q+1,1);
kernels{1} = vkern(0,'order',0,'domainincr',Ts);
for q = 1:Q
  k = zerokern(nLags,q);
  kernels{q+1} = vkern(k,'order',q,'domainincr',Ts);
end

ZeroKernels = kernels;



% 2. Generate kernels for each row and add them to the kernels



for i = 1:nparallel
  RowKernels = ZeroKernels;
  % initialize kernels with kernels of first 1 or 2 elements
  ss = subsys{i,1};
  if isa(ss,'kern')
    q = get(ss,'order');
    hlen = get(ss,'nLags');
    RowKernels{q+1} = pad_kernel(ss,nLags);
  elseif isa(ss,'irf')  
    k = zeros(nLags,1);
    hlen = get(ss,'nLags');
    k(1:hlen) = get(ss,'data');
    kobj = RowKernels{2};
    set(kobj,'data',k);
    RowKernels{2} = kobj;
  elseif isa(ss,'polynom');  
    q = get(ss,'order');
    set(ss,'type','power');
    mc = get(ss,'coef');
    for j = 0:q
      kobj = RowKernels{j+1};
      k = get(kobj,'data');
      k(1) = mc(j+1)/(Ts^j);
      set(kobj,'data',k);
      RowKernels{j+1} = kobj;
    end
  else 
    error('Unsupported element');
  end    
  % now deal with the rest of the elements in the row
  for j = 2:nseries
    ss = subsys{i,j};
    if isa(ss,'kern')
      RowKernels = kernel_convolve(RowKernels,ss);
    elseif isa(ss,'irf')
      RowKernels = kernel_convolve(RowKernels,ss);
    elseif isa(ss,'polynom');
      RowKernels = kernel_product(RowKernels,ss);
    else
      error('Unsupported element');
    end
  end
% Now add the kernels from this Row to the systems kernels.

  for q = 0:Q
    k = kernels{q+1};
    krow = RowKernels{q+1};
    set(k,'data', get(k,'data') + get(krow,'data'));
    kernels{q+1} = k;
  end
  
end % i loop

vs = vsin;
set(vs,'elements',kernels);



function kern = pad_kernel(k,nLags);
% pads kernel out to length nLags

hlen = get(k,'nLags');
if get(k,'nlags') < nLags
  kd = get(k,'data');
  q = get(k,'order');
  newkern =  zerokern(nLags,q);
  inds = 1:hlen;
  command = 'newkern(inds';
  for i = 2:q
    command = [command,',inds'];
  end
  command = [command,') = kd;'];
  eval(command);
  kern = k;
  set(kern,'data',newkern);
else
  kern = k;
end


function kernels = kernel_product(kerns,subsys)
% transforms kernels with a static nonlinearity.


% make sure the polynomial is a power series,
set(subsys,'type','power');
q = get(subsys,'order');
mc = get(subsys,'coef');

kernels = kerns;
%ks = get(kerns,'elements');
Q = length(kerns);


% Qmax will be the position of the last non-zero kernel
Qmax = 0;
for i = 1:Q
  ktest = double(kerns{i});
  if any(ktest(:))
    Qmax = i;
  end
end




new_ks = kerns;
for i = 1:Qmax
  new_ks{i} = 0*kerns{i};
end

% set the zero-order kernel equal to the zero order polynomial coefficient.
k0 = new_ks{1};
set(k0,'data',mc(1));
new_ks{1} = k0;

prod = kerns;
for i = 1:q
  % prod stores the kernels of the order i term.
  for j = 2:Q
    new_ks{j} = new_ks{j} + mc(i+1)*prod{j};
  end
  if i < q
    % generate the next order term
    prod = kproduct(prod,kerns);
  end
end

kernels = symmetrize_kernels(new_ks);
%set(kernels,'data',new_ks);

return


function kernels = kproduct(k1s,k2s);

% find the positions of the last non-zero kernels in k1s and k2s


% Qm1 will be the position of the last non-zero kernel in k1s
% Qm2 will be the position of the last non-zero kernel in k2s
Qm1 = 0;
Qm2 = 0;
Q = length(k1s);

for i = 1:Q
  ktest = double(k1s{i});
  if any(ktest(:))
    Qm1 = i;
  end
  ktest = double(k2s{i});
  if any(ktest(:))
    Qm2 = i;
  end  
end


 


kernels = k1s;
for i = 1:Qm1
  kernels{i} = 0*k1s{i};
end

for i = 2:Qm1
  k1 = double(k1s{i});
  k1 = k1(:);
  for j = 2:Qm2
    k = (i-1)+(j-1)+1;
    kern = 0*kernels{k};
    kk = double(kern);
    k2 = double(k2s{j});
    k2 = k2(:);
    kk = reshape(k1*k2',size(kk));
    set(kern,'data',kk);
    kernels{k} = kernels{k} + kern;
  end
end

%keyboard



function kernels = symmetrize_kernels(ks);

Q = length(ks);
for i = 1:Q
  k = ks{i};
  if i < 3
    % zero and first order kernels are already "symmetric"
    kernels{i} = k;
  elseif i == 3
    kdat = double(k);
    set(k,'data',(kdat+kdat')/2);
    kernels{i} = k;
  else
    % generate all possible permutations of i indeces
    % use PERMUTE to rearrange the kernel, and add'em up
    % divide by the number of permutations.
    perms = permutations(i-1);
    [np,i2] = size(perms);
    kdat = double(k);
    kk = 0*kdat;
    for j = 1:np
      kk = kk + permute(kdat,perms(j,:));
    end
    set(k,'data',kk/np);
    kernels{i} = k;
  end
end



function perms = permutations(i);

if i == 1 
  perms = 1;
else
  pmm1 = permutations(i-1);
  [nr,nc] = size(pmm1);
  perms = zeros(i*nr,nc+1);
  perms(1:nr,:) = [i*ones(nr,1), pmm1];
  for j = 2:i
    rowblock = [nr*(j-1)+1:nr*j]';
    perms(rowblock,:) = [pmm1(:,1:j-1) i*ones(nr,1) pmm1(:,j:nc)];
  end
end
    