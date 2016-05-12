% function D = select_model(M,pt)
% 
% This function operates on tvm objects
% This function allows the user to selct one model in the tvm object
% 
% pt = is the point to extract
% 
% This function return either an 'nlbl' or 'irf' object, depending on the
% original model type
% 
% V01-01    Nov 10/08   TSV  Initial code developed 
% V01-02    Nov 14/08   TSV  Extended to work on 'irf' model type 


function D = select_model(M,pt)

ModelType = get(M,'model_type');

switch ModelType
    case 'irf'
        
        temp  = get(M, 'data');  %get all models
        D = temp{pt}; 
         
    case 'nlbl'
               
        temp = get(M,'data');
        D = temp{pt};
        
    otherwise
        error (['tvm/select_model - model type:'  ModelType ' not defined']);
         
end
