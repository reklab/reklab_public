function P = reklabPaths (option);
%  P = reklabPaths (option, path_type
% set paths for reklab_public
% option: [add/remove]

% find base address to code can be used on any machine

if nargin < 1,
    option='Add';
end

s=(which('reklabPaths'));
k=findstr(s,[ filesep 'reklab_public' filesep]);
base=s(1:k);
disp(['Base path:' base]);

option=lower(option);

p= {'reklab_public/nlid_tools'  ...
    'reklab_public/nlid_tools/nlid_demo' ...
    'reklab_public/nlid_tools/nlid_util' ...
    'reklab_public/utility_tools' ...
    'reklab_public/stiffnessID' ...
    'reklab_public/smi-2.0-devel' };



for i = 1:length(p),
    p1=strcat(base,p{i});
    if strcmp (option,'add')
        disp(['Adding: ' p{i}]);
        addpath(p1);
    elseif strcmp(option,'remove')       
        disp(['Removing: ' p{i}]);
        rmpath(p1);      
    else
        error ([ 'Invald option for rkpaths:' option ]);
    end
    
end


