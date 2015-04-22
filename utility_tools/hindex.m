function h=hindex (oper, ht, varargin)
% hindex (oper, ht, varargin) haded index array operations
% $Revision: 1.11 $
% Usage:
% h=hindex (oper, ht, varargin)
%operations include:
% ht= hindex('init', n)  - initialize T as a hash table with initial size n
% i=hindex('put',ht, x) - put values into table; i=index values;
% i=hindex('get', ht,x) - get index for x
% x=hindex('val', ht, i) -return values for index i
% x=hindex('valuearray', ht ) - returns a cell  array of values by key. 
% 
% Note: the val operation is inefficient. UseValarray option for multiple
% calls. 
% 
% x may be a numeric value or array, a cell array, or a character string
%
%
% Example
%   H = hindex('init', 100)   % Initialize with room ofor 100 values
%   x= [ 100 200 300]
%   y= hindex('put', H, x) %
%   k =hindex('get',H, 200)  % Returns  k=2;
%   v =hindex('val',H, 3)    % return v = 300
% VA =hindex('valuearray',H  % VA = { 100 2000 300}
% $Revision: 1.11 $


switch lower(oper),
    case ('init');
        if nargin==2,
            h=java.util.HashMap(ht);
        else
            h=java.util.HashMap;
        end
    case 'put'
        keyVal=varargin{1};
        nVal=length(keyVal);
        if isnumeric(keyVal)
            argType='numeric';
        elseif iscell(keyVal),
            argType='cell';
        elseif ischar(keyVal),
            argType='char';
            nVal=1;
        else
            error ('Bad data type');
        end
        h=zeros(nVal,1);
        for i=1:nVal,
            switch argType
                case 'numeric'
                    curVal=keyVal(i);
                case 'cell'
                    curVal=keyVal{i};
                case 'char'
                    curVal=keyVal;
            end
            htemp=ht.get(curVal);
            if isempty(htemp),
                ht.put(curVal,ht.size+1);
                htemp=ht.size;
            end
                h(i)=htemp;
        end
    case 'get' % Return an index for a value
        keyVal=varargin{1};
        if ischar(keyVal),
            nVal=1;
        else
            nVal=length(keyVal);
        end
        h=zeros(nVal,1);
        for i=1:nVal,
            if ischar(keyVal),
                curVal=keyVal;
            elseif iscell(keyVal),
                curVal=keyVal{i};
            else
                curVal=keyVal(i);
            end
            temp=ht.get(curVal);
            if isempty(temp),
                h(i)=nan;
            else
                h(i)=temp;
            end
        end
    case 'val' % return a value for an index
        idx=varargin{1};
        idxList=ht.values;
        idxArray=idxList.toArray;
        idxArrayLen=length(idxArray);
        I=zeros(idxArrayLen,1);
        for i=1:idxArrayLen,
            I(i)=idxArray(i);
        end
        idxLen=length(idx);
        J=zeros(idxLen,1);
        for k=1:idxLen,
            J(k)=find(I==idx(k));
        end
        valList=ht.keySet;
        valArray=valList.toArray;
        h = cell(length(J),1);
        double_values = true;
        for j=1:length(J),
            h{j}=valArray(J(j));
            double_values = double_values & strcmp(class(h{j}),'double');
        end
        if(double_values)
            temp = h;
            h = zeros(size(temp));
            for j = 1:length(h)
                h(j) = temp{j};
            end
        end
        case 'valuearray' % Return an array of values by index
        idxArray=ht.values.toArray;
        valArray=ht.keySet.toArray;
        for i=1:ht.size,
            j=idxArray(i);
            h{j}=valArray(i);
        end
        numericFlag=all(cellfun(@isnumeric,h));
        if numericFlag,
            for i=1:length(h),
                temp(i)=h{i};
            end
            h=temp;
        end
            

    otherwise
        error(['Operation not defined:' oper ]);
end