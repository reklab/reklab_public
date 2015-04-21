function [W, flags] = multi_test(V, max_order);
%
%  usage W = multi_herm2(V, max_order);
%
%  given a collection of column vectors, V, this function returns
%  a matrix of all of the Tcheb functions up to max_order applied
%  to all of the vectors in V.

% Copyright 1999-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

if max_order==0,
   W=V*0 +1;
   return;
end




[nr,nc] = size(V);


%  create a matrix of flags, such that flags (i,j) points to the 
%  column of the i'th basis function raised to the j'th power.

flags = zeros(nc,max_order);
flags(:,1) = [2:nc+1]';
for order = 2:max_order
   % first column is offsetr by 1 from last value of previous order
   flags(1,order) = flags(nc,order-1)+1;
    for i = 2 : nc
        num_terms = flags(nc,order-1)-flags(i-1,order-1)+1;
        flags(i,order) = flags(i-1,order)+ num_terms;
      end
   end

   

% pre-allocate the W matrix
W = zeros(nr,flags(nc,max_order));
%  generate the single basis function Hermite polynomials and place 
%  them  in the correct columns of W

Te = zeros(nr,max_order+1);
Te(:,1) = ones(nr,1);

%  generate all of the functions involving a Hermite polynomial
%  applied to a single basis vector

for i = 1:nc
    Te(:,2) = V(:,i);
    for j = 2:max_order
        Te(:,j+1) = 2*V(:,i).*Te(:,j) - Te(:,j-1);
      end
    for j = 1:max_order
        W(:,flags(i,j)) = Te(:,j+1);
      end
  end
W(:,1) = ones(nr,1);

clear Te V


%  generate all of the functions involving a powers of a single input vector


  
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
                disp([index flags(v1,v1_order) j])
              end
          end
      end
  end




