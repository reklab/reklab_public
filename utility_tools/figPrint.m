function figPrint(figNum, fileName, printFlag, newFileFlag, pageOrient)
%  figPrint(figList, fileName)
% utility routine to contorl printinig of figures
% printFlag - logical variable; prnt if true
% figNum - figure to print
% fileName - name of pring file
% newfileFlag - create new file if true, otherwise append
if nargin<4,
    newFileFlag=false;
    
end
if nargin <5,
    pageOrient = 'portrait';
end
if nargin,3,
    printFlag=true;
end

if ~printFlag,
    return
end

if newFileFlag
    orient(figNum,pageOrient); 
    print (figNum, fileName, '-dpsc' ,'-bestfit');
    
else
    print (figNum, fileName, '-dpsc','-append','-bestfit');
    
end
