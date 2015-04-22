function S = subset_fun ( OP, SUBSET, setid, elements )
% subset_fun - supports subset operatins
% $Revision: 1.2 $
% S = subset_fun  ( OP,  SUBSET, setid, elements ) $Revision: 1.2 $
% OP - operations
%       "ADD" - add a subset
%       "ADD_ELEMENTS" - add elements to a set
%       "ADJACENT" - return adjacency matrix for defined subets
%       "FIND_SUPERSET" - find ids containing all elements.
%       "ELEMENTS" - return all elements in sets
%       "INDEX" - returns index to subset elemets
%       "HELP"
%       "INIT" - initialze set
%       "INTERSECT - interection of two sets
%       "MEMBERS" - returns member of elements insetid
%       "SETDIFF" - difference between two sets
%       SUBSET_ELEMENTS - return memers of  subset setid
%       SUBSET_INDEX - return SUSET index for subsetid
%       "SET" - return set matrix
%       "UNIQUE" - returns elements of setid which are unique to the subset

%
% u = set element(s);
% v set element
% fucntions to support subsets
persistent SUBSET_ELEMENTS
if nargin==0,
    subset_fun('help');
    return
end

switch upper (OP),
    case 'ADD' % Add a subset
        i=length([SUBSET.ID])+1;
        SUBSET(i).ID=setid;
        SUBSET(i).ELEMENTS=elements;
        SUBSET_ELEMENTS= union (SUBSET_ELEMENTS, elements);
        S=SUBSET;
    case 'ADD_ELEMENTS' % Add a elements to a subset
        i=find([SUBSET.ID]==setid);
        if length(i)==0,
            S=subset_fun('add', SUBSET, setid,elements);
            return
        elseif length(i)==1,

            SUBSET(i).ELEMENTS=union(SUBSET(i).ELEMENTS, elements);
            SUBSET_ELEMENTS= union (SUBSET_ELEMENTS, elements);
            S=SUBSET;
        else
            error ('subset index multiply defined');
        end
    case 'ADJACENT' % Return adjacency matix for subsets
        nset=length(SUBSET);
        for i=1:nset,
            ei=SUBSET(i).ELEMENTS;
            for j=1:nset,
                ej=SUBSET(j).ELEMENTS;
                S(i,j)=length(intersect(ei,ej));
            end
        end
    case 'FIND_SUPERSET' % find subset contains all elements
        S=[];
        for i=1:length(SUBSET),
            if isempty(setdiff(elements,[SUBSET(i).ELEMENTS])),
                S=[ S SUBSET(i).ID ];
            end
        end
            
    case 'INDEX' % Return index to subset elements
        S=[];
        for i=1:length(SUBSET(setid).ELEMENTS),
            S=[ S find(SUBSET(setid).ELEMENTS(i)== SUBSET_ELEMENTS)];
        end


    case 'HELP'
        hlp = {  ' R = subset_fun ( OP, u, v  )' ...
            'defined operations:' ...
            'COLLECTION -  Return the collection of sets'  ...
            };
        for i=1:length(hlp),
            disp(hlp{i});
        end
    case 'INIT' % Initialize set collection
        clear SUBSET;
        SUBSET.ID = [];
        SUBSET.ELEMENTS=[];
        SUBSET_ELEMENTS=[];
        S=SUBSET;
    case 'INTERSECT'
        S=intersect(SUBSET(setid).ELEMENTS, SUBSET(elements).ELEMENTS);
    case 'ELEMENTS' % Return all set elements
        S=SUBSET_ELEMENTS;
    case 'MEMBERS' % Return members of "elements" in setid
        i=find(ismember(elements, SUBSET(setid).ELEMENTS));
        S=elements(i);
    case 'SET'
        S=SUBSET;
    case 'SETDIFF'
        S=setdiff(SUBSET(setid).ELEMENTS, SUBSET(elements).ELEMENTS);
    case 'SUBSET_ELEMENTS'
          i=find([SUBSET.ID]==setid);
          S=SUBSET(i).ELEMENTS;
    case 'SUBSET_INDEX' % REturns index nuymber of susbset setId;
        S=find([SUBSET.ID]==setid);
    case 'UNIQUE'
        iset = find([SUBSET.ID]==setid);
        if length(iset)==0,
            S=[];
            return
        end
        S=SUBSET(iset).ELEMENTS;
        for i=1:length(SUBSET),
            if SUBSET(i).ID ~= setid,
                S=setdiff(S,SUBSET(i).ELEMENTS);
            end
        end
    otherwise
        error (['Undefined operation: ' OP]);
end

return
