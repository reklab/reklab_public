function N = subsref (N,S);
% overlaid subsref for tvm objects
% N{i} returns all models elements at one time index i 
% N(i,j) returns i,jth model element at all times 

for icount=1:length(S),
    
if strcmp (S(icount).type, '.'),
    switch lower(S(icount).subs)
    case 'data'
        N=get(N,'data');
    otherwise
        disp('not defined');
    end
else
    if iscell(N),
        N=subsref(N,S(icount:end));
    else
    
    
    
    E=get(N,'data');
    [nx(1),nx(2)]=size(E);
    for i=1:2,
        ix{i}=1:nx(i);
    end
    for i = 1:length(S(icount).subs),
        if ischar(S(icount).subs{i}),
            ix{i}=1:nx(i);
        else
            ix{i}=S(icount).subs{i};
        end
    end
    
    if strcmp (S(icount).type, '{}'),
        s=ix{1};
        d=domain(N);
        for j=1:length(s),
            ed{j,1}=E{s(j)};
            dnew(j)=d(s(j));
        end
        set(N,'data',ed,'domainvalues',dnew);
        
        
        % () return i,jth element of model
    elseif strcmp (S(icount).type, '()'),
        for j=1:length(E),
            et=E{j};
            ed{j,1}=et{ix{1},ix{2}};
        end
        
        set (N,'data',ed,'Model_Type',class(ed{1,1}));
        
        
    else
        N=NaN;
        error ('NLM subsref not yet implemented for this syntax');
    end
end
end
end

return
% tvm/subsref
