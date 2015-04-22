function figPrint(figList, fileName)
%  figPrint(figList, fileName)
% print multiple figures to a file 
if nargin==0,
   figLIst=get(0,'children');
end
for i=1:length(figList),
    figure(figList(i));
   if i==1,
           eval( ['print -dpsc '  fileName ';' ] );
        else
           eval( [' print -dpsc -append ' fileName ';' ]);
        end
end
end
