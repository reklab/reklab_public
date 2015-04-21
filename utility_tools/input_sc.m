function ans = input_sc ( prompt,default,options)
% input a string with optional  default and list
% 	function val = input_s ( prompt,default,options)
%
%       prompt - prompt string
%       default - default value
%       optons  - cell array list of acceptable options
%

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

% 28 Aug 98 Problem with "\" chars in promt
% 2001 Oct Modify to use cells 
if nargin ==1,
  default = '?';
else
  x = default;
end
done=0;
while (~done),
   ptemp = [ prompt '  [' default  ']:'];
   ptemp=strrep(ptemp,'\','\\');
  temp = input (ptemp,'s');
  if isempty (temp),
    ans = default;
  else 
    ans=temp;
  end
  if (nargin==3),
    done = strmatch(ans, options);
    if (length(done)==0),
      disp('Invalid. Options are:');
      for i=1:length(options),
          disp(options{i})
      end;
    end
  else
    done =1;
  end
end

