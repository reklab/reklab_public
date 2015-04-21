function ans = input_l ( prompt,default)
% input logical y/n answer
% 	function val = input_d ( prompt,default)

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

defstr='N';
if nargin ==1,
	default = 'y';
else
	if(default),
		defstr='Y';
	end

end
temp = input ([ prompt ' (y/n) [' defstr  ']:'],'s');                         
if isempty (temp),
	val = defstr;
else 
	val=temp;
end
if ( val(1) == 'y') | (val(1)=='Y'),
	ans = 1;
else
	ans = 0;
end

