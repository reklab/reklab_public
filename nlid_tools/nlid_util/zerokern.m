function kern = zerokern(hlen,q);
% ZEROKERN - generates an order q array full of zeros.
%
% $Revision: 1.4 $

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


switch q
  case 0
    kern = 0;
  case 1
    kern = zeros(hlen,1);
  case 2
    kern = zeros(hlen,hlen);
  otherwise
    kern = zeros(hlen^q,1);
    % reshape the array into a square, q dimensional tensor.
    % Since it must work for all values of q, we have to generate the command
    % with a loop, and then use eval to execute it.
    command = 'kern = reshape(kern,hlen';
    for order = 2:q
      command = [command,',hlen'];
    end
    command = [command,');'];
    eval(command);
end

