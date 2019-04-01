function h = dockFig(mode, figNum)
% dockFig - dock/undock figure(1)
% mode - dock/undock
% figNum - list of figures to dock/undock. If not specifed apply to all
% figures
if nargin==1,
    figNum=get(0,'children');
end
for i=1:length(figNum),
    switch mode
        case 'dock'
            set(figNum(i),'windowstyle','docked'); 
        case 'undock'
             set(figNum(i),'windowstyle','normal'); 
        otherwise
            error(['Bad mode command:' mode]);
    end
end


