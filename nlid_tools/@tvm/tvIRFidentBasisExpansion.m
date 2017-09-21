

end
function [H, x_pred, Extra] = np_TV_ident(X, Y, Basis, varargin)

% Identifies dynamics of linear, time-varying, non-parametric system.
% Diegos method
% [H, x_pred] = np_TV_ident(X, Y, varargin)
%
% X:              input matrix
% Y:              output matrix
% Basis:          Basis functions
% Additional arguments:
% 'nLags'         Length of linear component IRF in samples
% 'nSides'        1 sided or 2 sided IRF
% 'domainIncr'    Sampling rate in seconds
% 'peridic'       Is the data periodic {'yes' / 'no'}
% 'method'        Linear identification method {'Bayes' / 'OLS'}
%
% H:              matrix of identified TV-impulse response function
% x_pred:         predicted output. Free run of the estimated model
%
% The matrices X and H have the following format:
%
% X:  ==== Trials ====>      H:  ==== Lag Time ====>
%          ||                         ||
%          ||                         ||
%          ||                         ||
%       Discrete                   Discrete
%         Time                       Time
%          ||                         ||
%          ||                         ||
%          ||                         ||
%          ||                         ||
%          \/                         \/
%
% The format of the output matrix Y is identical to that of X.
%




[Ns, trials]= size(X);

if length(Y(:,1))~= Ns || length(Y(1,:))~= trials
    error('Input and output matrices should be of equal size')
end

if length(Basis)~=Ns
    error('The basis functions should be of equal lenght as the data')
end

if nSides==1
    posLags=nLags;
    negLags=0;
elseif nSides==2
    posLags=nLags;
    negLags=-nLags;
else
    error('nSides should be equal to 1 or 2')
end

dt=domainIncr;


%creating the regressor matrix
%The type of matrix depends if the data is periodic or not and if the IRF
%is one sided or two sided
BASIS_USED=[];
n_bas=size(Basis,2)-1;

if strcmp(periodic,'yes')
    for p=1:trials
        
        aux=1;
        for i=negLags:0
            if p==trials
                %REG_M{p}(:,aux)=[pos(-i+1:end,p);zeros(-i,1)];
                REG_M{p}(:,aux)=[X(-i+1:end,p);X(1:-i,1)];
                aux=aux+1;
            else
                REG_M{p}(:,aux)=[X(-i+1:end,p);X(1:-i,p+1)];
                aux=aux+1;
            end
            
        end
        for i=1:posLags
            if p==1
                %REG_M{p}(:,aux)=[zeros(i,1);pos(1:end-i,p)];
                REG_M{p}(:,aux)=[X(end-i+1:end,end);X(1:end-i,p)];
                aux=aux+1;
            else
                REG_M{p}(:,aux)=[X(end-i+1:end,p-1);X(1:end-i,p)];
                aux=aux+1;
            end
            
        end
        BAS{p}=[];
        for k=1:size(REG_M{p},2)
            BAS{p}=[BAS{p} Basis.*repmat(REG_M{p}(:,k),1,n_bas+1)];
        end
        %staking all the matrices together
        BASIS_USED=[BASIS_USED;BAS{p}];
    end
elseif strcmp(periodic,'no')
    
    for p=1:trials
        
        aux=1;
        for i=negLags:0
            if p==trials
                REG_M{p}(:,aux)=[X(-i+1:end,p);zeros(-i,1)];
            end
            
        end
        for i=1:posLags
            if p==1
                REG_M{p}(:,aux)=[zeros(i,1);X(1:end-i,p)];
            end
            
        end
        BAS{p}=[];
        for k=1:size(REG_M{p},2)
            BAS{p}=[BAS{p} Basis.*repmat(REG_M{p}(:,k),1,n_bas+1)];
        end
        %staking all the matrices together
        BASIS_USED=[BASIS_USED;BAS{p}];
    end
else
    error('''periodic'' can only be ''yes'' or ''no''')
end


%estimation of TV-IRF
if strcmp(method,'Bayes')
    [~,mu,noise_var,lk,COV]=TV_Bayes(vec(Y),[],BASIS_USED);
    
elseif strcmp(method,'OLS')
    mu=pinv(BASIS_USED)*vec(Y);
    noise_var=var(vec(Y)-BASIS_USED*mu);
    TEMP=BASIS_USED'*BASIS_USED;
    COV=(noise_var/Ns)*(eye(size(TEMP))/(TEMP));
    lk=NaN;
    
else
    error('Only ''Bayes'' and ''OLS'' ID methods are implemented')
end

%function output
x_pred=zeros(size(X));
H=zeros(Ns,size(REG_M{1},2));
for p=1:trials
    x_pred(:,p)=BAS{p}*mu;
end
aux=1;
for k=1:size(REG_M{1},2)
    H(:,k)=(1/dt)*Basis*mu(aux:aux+n_bas);
    aux=aux+n_bas+1;
end

Extra.id_params=mu;
Extra.Covariance_id_params=COV;
Extra.noise_variance=noise_var;
Extra.log_likelihood=lk;
return

end

