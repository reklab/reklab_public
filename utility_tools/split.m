function SP = split ( S, delimit)
% split a string into components based on delimiter
% $Revision: 1.1 $
% Two consecutive delimiters returns a Nab
%Ignore leading delimiters
% Return original string if delimiter is not found; 
idelimit=strfind(S,delimit);
if length(idelimit)==0,
    SP{1}=S;
    return
end
istart=1;
j=1;
for i=1:length(idelimit)
    istop=idelimit(i); 
    if istop==1,
        istart=length(delimit)+1;
    elseif istop==istart,
        istart=istop+1;
        SP{j}='';
        j=j+1;
    else
        SP{j}=S(istart:istop-1);
        istart=idelimit(i)+length(delimit);
        j=j+1;
    end
end
if length(S)>=istart,
    SP{j}=S(istart:end);
end

return


