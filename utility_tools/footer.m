function handle = footer(string)
% printer a footer message at the bottom of the screen
axes('position',[0 0 1 1],'Box','off','Visible','off');
handle=text (0.01, 0.02,string, 'units','normalized','fontsize',8, ...
    'fontangle','italic');

