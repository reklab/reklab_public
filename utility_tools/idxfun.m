function  JDX = idxfun ( op, IDX, x );
%  idxfun -  support indexed arrays for clusters
% $Revision $
% JDX = idxfun ( op, IDX, x );
% op -  "init" - initialize arrays
%       "add" add x to index
%       "index" return index values for x
%       'invert - return the inverted index
% IDX - indexed array
% x - value on which to operate
% JDX(x,1)- contains unique index for x
% JDX(x,2) - contains number of times value has been called
xdimMax = 10e7;
switch  lower(op),
    case 'init'
        JDX=spalloc(xdimMax,2,10000);
    case 'add'
        for i=1:length(x),
            ix=x(i);
            if (ix > xdimMax)
                error (['Index value of ' num2str(ix) ... 
                        ' exceeds maximum allowable value of ' num2str(xdimMax)]);
            end
            if IDX(ix,1)==0,
                IDX(ix,1)=max(IDX(:,1))+1;
            end
            IDX(ix,2)=IDX(ix,2)+1;
        end
        JDX=IDX;
      case 'invert' % Returns index values of x
        i=find(IDX(:,1)>0);
        for j=1:length(i),
            JDX(j,1)=i(j);
            JDX(j,2)=IDX(i(j),2);
        end
    otherwise
        error(['Invalid option: ' op]);
end


