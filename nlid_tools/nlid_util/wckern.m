function kernel = wckern (g,Q);
% WCKERN - generates kernel of order Q for wiener system with g as its IRF.
%
% $Revision: 1.3 $

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 



switch Q
  case 0
    kernel = 1;
  case 1
    kernel = g;
  case 2
    g = g(:);
    kernel = g*g';
  otherwise
    % for kernels of order 3 or greater, we will need to use a multidimensinal
    % array.  First, create an empty kernel of the appropriate dimensions.     

    hlen = length(g);
    indeces = hlen*ones(1,Q);
    kernel = zeros(indeces);


    % the idea here, is squish the existing kernel into a column, multiply by
    % the transpose of the IRF, then squish the result into a column...
    % do this until the column has length hlen^Q

    g = g(:);
    gt = g';
    k = g*gt;
    k = k(:);
    for order = 3:Q
      k = k*gt;
      k = k(:);
    end

    % reshape the result into a square, Q dimensional tensor.
    % Since it must work for all values of Q, we have to generate the command
    % with a loop, and then use eval to execute it.
    command = 'kernel = reshape(k,hlen';
    for order = 2:Q
      command = [command,',hlen'];
    end
    command = [command,');'];
    eval(command);
end

