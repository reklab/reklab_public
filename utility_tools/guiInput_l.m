function ans = guiInput_l ( prompt,default)
% input logical y/n answer
% 	function val = input_d ( prompt,default)

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

defstr='No';
if nargin ==1,
	default = 'Yes';
else
	if(default),
		defstr='Yes';
    else
        defstr='No';
    end

end
drawnow
val=questdlg(prompt, 'Logical input', defstr);
drawnow
switch val
    case 'Yes'
        ans=true;
    case 'No'
        ans=false;
    otherwise
            error('Input cancelled');
end
        

