% function D = select_nlbl(hammer,pt)
% 
% This function allows the user to selct one NLBL element in the tvm object
% This function returns an NLBL object
% pt = is the time point to extract
% 
% V01-01    Nov 10/08   TSV  Initial code developed 


function D = select_nlbl(hammer,pt)

D = nlbl; %create an NLBL object
temp = get(hammer,'data');
set(D,'elements', temp(pt));

end
