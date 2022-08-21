function  printFigures(figureList, fileName)
% printFigures (figureList, fileName)


nFig=length(figureList)
appendFlag=false;
for i=1:nFig
    figure(figureList(i));
    exportgraphics(gcf,fileName,"Append",appendFlag)
    appendFlag=true;
end