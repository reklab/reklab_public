% function C = select_irf(hammer,pt)
% 
% This function allows the user to extract one IRF from an NLBL element in 
% the tvm object
% This function returns an IRF object
% pt = is the time point to extract
% 
% V01-01    Nov 10/08   TSV  Initial code developed 


function C = select_irf(hammer,pt)

A = get(hammer, 'data');
B = get(A{pt},'elements');
C = B{1,2};

end
