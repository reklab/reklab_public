% function hammer_short = extract_nlbl(hammer, start, final)
% 
% This function removes points from the beginning and the end of the tvm
% object of model type 'nlbl'
% start = number of nlbl elements to remove at the beginning
% final = number of nlbl elements to remove at the end
% 
% V01-01    Nov 10/08   TSV  Initial code developed 


function hammer_short = extract_nlbl(hammer, start, final)

hammer_short = hammer;  %create identical object
D = get(hammer,'data');
D = D(start+1:end-final);   %retain wanted points
set(hammer_short,'data',D); 

end