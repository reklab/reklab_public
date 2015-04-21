function val = guiInput_d ( prompt,default,xmin,xmax)
% guiInput_d  gui based interqactive input of numeric values with defaults and  range checking
% 	function val = input_d ( prompt,default,[xmin],[xmax])
% 
%  Users may enter a numeric value or a matal expression which will be
%  evaluated in the calling workspace
% $Revision: 1.5 $
% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if (nargin < 4),
	xmax = 1.e20;
end
if (nargin < 3),
	xmin = -1e20;
end
if nargin < 2,
	x = 0;
else
	x = min(default, xmax);
end
done = false;
while (~ done),
    defaultValue= sprintf(' %d',x);
    drawnow;
	temp = inputdlg ( char(prompt) ,'Number', 1,{defaultValue});
    drawnow
	if isempty (temp),
		error('Input cancelled');
    else
       temp=evalin ('base',temp{1});
		if any(temp > xmax),
			disp(['Too big. Values must be < ' num2str(xmax)])
		elseif any(temp < xmin),
			disp(['Too small. Values must > ' num2str(xmin)])
		else
			val = temp;
			done = true;
		end
	end
end
