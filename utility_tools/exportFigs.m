function exportFigs (figList, exportName, fFormat, resolution )
% export one or more figures to a file
% figLIst - list of figures to export
if nargin<4
    resolution=300;
end
if nargin<3
    fFormat='JPG';
end
if nargin<2
    exportName='Figure';
end
nFig=length(figList);
for i=1:nFig
    curFig=figList(i);
    f=figure(curFig)
    curFileName= [ exportName '_' num2str(curFig) '.'  fFormat];
    exportgraphics(f,curFileName,'Resolution',resolution);
    disp(['Figure ' num2str(curFig) ' exported as ' curFileName]);
end

