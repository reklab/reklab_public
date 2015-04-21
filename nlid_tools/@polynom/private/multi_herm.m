function [W, flags] = multi_herm2(V, max_order);
%
%  usage W = multi_herm2(V, max_order);
%
%  given a collection of column vectors, V, this function returns
%  a matrix of all of the hermite functions up to max_order applied
%  to all of the vectors in V.

% Copyright 1998-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

[nr,nc] = size(V);

%  create a matrix of flags, such that flags (i,j) points to the 
%  column of the i'th basis function raised to the j'th power.

flags = zeros(nc,max_order);
flags(:,1) = [2:nc+1]';
for order = 2:max_order
  flags(1,order) = flags(nc,order-1)+1;
  for i = 2 : nc
    num_terms = flags(nc,order-1)-flags(i-1,order-1)+1;
    flags(i,order) = flags(i-1,order)+ num_terms;
  end
end


% pre-allocate the W matrix
if max_order > 0
  W = zeros(nr,flags(nc,max_order));
else
  W = zeros(nr,1);
end

%  generate the single basis function Hermite polynomials and place 
%  them  in the correct columns of W


He = zeros(nr,max_order+1);
He(:,1) = ones(nr,1);

%  generate all of the functions involving a Hermite polynomial
%  applied to a single basis vector

for i = 1:nc
  He(:,2) = V(:,i);
  for j = 2:max_order
    He(:,j+1) = V(:,i).*He(:,j) - (j-1)*He(:,j-1);
  end
  for j = 1:max_order
    W(:,flags(i,j)) = He(:,j+1);
  end
end
W(:,1) = ones(nr,1);

clear He V

%  Now, using the functions that we just created, fill in the 
%  rest of the matrix



for order = 2:max_order
  for v1 = 1 : nc-1
    index = flags(v1,order);
    for v1_order = order-1:-1:1
      term1 = W(:,flags(v1,v1_order));
      rem_order = order - v1_order;
      
      %           find the terms of order rem_order, whose leading variable
      %           is 'greater' than v1

      first_term = flags(v1+1,rem_order);
      last_term = flags(nc,rem_order);
      for j = first_term:last_term
	index = index+1;
	W(:,index)=term1.*W(:,j);
	%                disp([index flags(v1,v1_order) j])
      end
    end
  end
end




