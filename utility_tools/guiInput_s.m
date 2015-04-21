function ans = guiInput_s ( prompt,default,options)
% input a string from guiwith optional  default and list
% 	function val = input_s ( prompt,default,options)
%
%       prompt - prompt string
%       default - default value
%       optons  - list of acceptable options
%

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

% 28 Aug 98 Problem with "\" chars in promt
% 2001 Oct Modify to use cells 
if nargin ==1,
  default = '';
else
  x = default;
end
done=0;
while (~done),

   ptemp=strrep(prompt,'\','\\');
   drawnow
  temp = inputdlg(char(ptemp),'Enter a string',[1 80],{default});
  drawnow
  if isempty (temp),
    ans = error('input cancelled');
  else 
    ans=temp{1};
  end
  if (nargin==3),
    done = ismember(options,ans);
    if ~done,
      disp('Invalid. Options are:');
      disp(options);
    end
  else
    done =1;
  end
end

