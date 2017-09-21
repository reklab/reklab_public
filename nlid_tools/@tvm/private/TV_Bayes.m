function [x_est,mean_pars,noise_var,lk,COV,par,x_est1]=TV_Bayes(y,u,BAS)

%This algorithm implements the Relevance Vector Machine for selection of
%the best basis function for Time-Varying identification. 

%For reference see:
%Michael E. Tipping,'Sparse Bayesian Learning and the Relevance Vector
%Machine',Journal of Machine Learning Research 1 (2001) 211{244

if ~isempty(u)
    [N,nb]=size(u);
    BASIS=[];

    if nargin<3
        for i=1:nb
        BAS{i}=ones(N,1);
        end
    end
    for i=1:nb
        N_basis(i)=size(BAS{i},2);
        BASIS=[BASIS BAS{i}.*repmat(u(:,i),1,N_basis(i))];
    end
else
    BASIS=BAS;
    N=length(y);
end

par=BASIS\y;
%par=robustfit(BASIS,y,[],[],'off');
x_est=BASIS*par;
x_est1=x_est;
error=y-x_est;
beta=1/var(error);

alpha=(1/beta)*diag(eye(size(BASIS,2))/(BASIS'*BASIS));
alpha=length(par)/(2*(par'*par));

mean_pars=par';

Max_iter=50;
clear lk


alpha=1/(10^-10);
mean_pars=randn(size(mean_pars));

for iter=1:Max_iter

    Be=beta*(BASIS'*BASIS);
    lambda=eig(Be);

    gamma=sum(lambda./(alpha+lambda));
    
    alpha=gamma/(mean_pars*mean_pars');
    
    error=y-BASIS*mean_pars';
    
    beta=1/((1/(N-gamma))*(sum(error.^2)));
    
    A=beta*(BASIS'*BASIS)+alpha*eye(size(BASIS,2));
    
    mean_pars=(beta*(A\(BASIS'*y)))';
    
   lk(iter)=(size(BASIS,2)/2)*log(alpha)+(size(BASIS,1)/2)*log(beta)-...
    (beta/2)*norm(y-BASIS*mean_pars')^2-(alpha/2)*(mean_pars*mean_pars')-...
    (1/2)*log(det(A))-(N/2)*log(2*pi);

   if iter>1
        if((abs(lk(iter)-lk(iter-1)))/abs(lk(iter-1)))<1e-7
            break
        end
    end
    
end
mean_pars=mean_pars';
x_est=BASIS*mean_pars;
noise_var=1/beta;
COV=(eye(size(A))/(A));

return
