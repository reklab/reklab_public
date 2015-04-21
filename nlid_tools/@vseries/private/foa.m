function vs = foa(vs,u,y,Ts)

% primary entry point for FOA computation of Volterra kernels 
% 
% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 
hlen = get(vs,'nLags');
Qmax = get(vs,'vsOrderMax');
kernels = get(vs,'elements');
vk0 = kernels{1};
set(vk0,'nLags',hlen,'domainIncr',Ts);
N = length(u);
M = prod([hlen+1:hlen+Qmax])/factorial(Qmax);
%  M is the total number of regressors
[PmPr,G] = foa_init(u,y,hlen,Qmax);
[D,ord] = chol(PmPr);
if ord > 0
  disp('Hessian is poorly conditioned -- using SVD');
  Hinv = pinv(PmPr);
  A = Hinv*G;
else
  A = D\((D'\G));
end
%  and of course, the kernel values.
kernels = cell(Qmax+1,1);
h0 = A(1);
kernels{1} = vkern(vk0,'kernOrder',0,'dataSet',h0);
if Qmax > 0
  h1 = A(2:hlen+1);
  kernels{2} = vkern(vk0,'kernOrder',1,'dataSet',h1/Ts);
end
if Qmax > 1
  m = hlen + 1;
  for i = 1:hlen
    for j = i:hlen
      m = m + 1;
      if i == j
	h2(i,i) = A(m);
      else
	h2(i,j) = A(m)/2;
	h2(j,i) = A(m)/2;
      end
    end
  end
  kernels{3} = vkern(vk0,'kernOrder',2,'dataSet',h2/(Ts^2));
end
set(vs,'elements', kernels);
set(vs,'comment','Identified using Fast Orthogonal Algorithm');

% and compute the kernel variances
if 0
if nargout > 3
  if ord == 0
    Di = inv(D);
    Avar = (1/N)*diag(Di*Di');
  else
    Avar = (1/N)*diag(Hinv);
  end

  

  

  h0var = Avar(1);

  h1var = Avar(2:hlen+1);



  m = hlen + 1;

  for i = 1:hlen

    for j = i:hlen

      m = m + 1;

      if i == j

        h2var(i,i) = Avar(m);

      else

        h2var(i,j) = Avar(m)/4;

        h2var(j,i) = Avar(m)/4;

      end

    end

  end

  ssresid = sum(y.^2)-N*(A'*PmPr*A);

  gain = ssresid/(N-M);

  h0var = gain*h0var;

  h1var = gain*h1var;

  h2var = gain*h2var;

end

end



















