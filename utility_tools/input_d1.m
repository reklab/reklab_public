function val = input_d1 ( handle,default,xmin,xmax)
% decode a  string with default values and option range checking
%

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
	x = default;
end
done = 0;
	temp = str2num(get(handle,'string'));                         
	if isempty (temp),
		val = x;
		done = 1;
	else
		if (temp > xmax),
			disp('input_d1: Value to large. Setting to max');
			set(handle,'String',num2str(xmax));
			val =xmax;
		elseif (temp < xmin),
		  	disp('input_d1: Value too small. Setting to min');
			set(handle,'String',num2str(xmin));
			val =xmin;
		else
			val = temp;
			done = 1;
		end
	end

