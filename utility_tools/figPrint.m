function figPrint(figNum, fileName, printFlag, newFileFlag, pageOrient, formatType)
%  figPrint(figNum, fileName, printFlag, newFileFlag, pageOrient, formatType)
% utility routine to control printinig of figures
% printFlag - logical variable; prnt if true
% figNum - figure to print
% fileName - name of pring file
% newfileFlag - create new file if true, otherwise append
% formatType - output format as defined in matlab documentation [-dpsc ]
if nargin<4,
    newFileFlag=false; 
end
if nargin <6,
    formatType = '-dpsc';
end
if nargin <5,
    pageOrient = 'portrait';
end
if nargin<3,
    printFlag=true;
end

if ~printFlag,
    return
end

if newFileFlag
    orient(figNum,pageOrient); 
    print (figNum, fileName, formatType);
    
else
    print (figNum, fileName, formatType,'-append','-bestfit');
    
end
