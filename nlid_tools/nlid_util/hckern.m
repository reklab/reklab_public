function kernel = hckern(g,Q);
% HCKERN - generates Q'th order Volterra kernel of a Hammerstein cascade.
% $Revision: 1.3 $
% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

switch Q
  case 0
    kernel = sum(g);
  case 1
    kernel = g;
  case 2
    kernel = diag(g);
  otherwise
    % for kernels of order 3 or greater, we will need to use a multidimensinal
    % array.  First, create an empty kernel of the appropriate dimensions.     

    hlen = length(g);
    indeces = hlen*ones(1,Q);
    kernel = zeros(indeces);

    % There should be a much more efficient way of doing this,
    command = 'kernel(i,i';
    idx = [',i'];
    for j = 3:Q
      command = [command, idx];
    end
    command = [command,') = g(i);'];
    for i = 1:hlen
      eval(command);
    end
end

