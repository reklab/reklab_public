% function short = extract(M, start, final)
% 
% This function operates on tvm object
%   model type: 'nlbl' and 'irf'
% This function removes models at time points from the beginning and the 
% end of the tvm object 
% 
% start = number of  models to remove at the beginning
% final = number of models to remove at the end
% 
% V01-01    Nov 14/08   TSV  Initial code developed 


function short = extract(M, start, final)

short = M;  %create identical object
D = get(M,'data');  %extract all the models
D = D(start+1:end-final);   %retain wanted models
set(short,'data',D); 

end