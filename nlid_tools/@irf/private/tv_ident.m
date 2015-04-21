function [Hident,bound,sing_vectors,cpu] = tv_ident(X,Y,dt,arg1,arg2,approach,conf_level)
 
% Identifies dynamics of linear time-varying system.  
%
% [Hident,bound,sing_vectors,fl,cpu] = tv_ident(X,Y,dt,arg1,arg2,approach,conf_level)
%
% X:              input matrix 
% Y:              output matrix
% dt:             sampling interval
% arg1 and arg2:  arg1 = - length of (i.e. number of points in) anticipatory
%                        component of IRF
%                 arg2 = length of (i.e. number of points in) memory 
%                        component of IRF
%                                    OR
%                 arg1 = 'one' (1-sided IRF) or 'two' (2-sided IRF)
%                 arg2 = length of (i.e. number of points in) the IRF
% approach:       'tvfil' - finds least-squares solution using the data itself.
%                 'corr' - finds least-squares solution using correlation functions.
%                 'pseudo' - correlation function approach as with 'corr', but
%                            uses Dave's technique (adapted to tv case) to select singular vectors
%                 Default is 'corr'.
% conf_level:     confidence level for IRF bound (only applies when approach is 'pseudo')
%                 If conf_level is NaN, the confidence bounds are not computed.  Default is NaN.
%
% Hident:         matrix of identified dynamics
% bound:          matrix of confidence bounds on identified dynamics
% sing_vectors:   number of singular vectors used to compute the inverse or pseudoinverse
% cpu:            cputime required to perform identification
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

% Copyright 1998-2003, Robert E Kearney, David T Westwick and Mireille Lortie
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

% Check that number of input arguments is valid.
num_args = nargin;
error(nargchk(5,7,num_args));

% Check that input arguments are valid.
if ~(isnumeric(X))
   error('The first argument must be the input matrix.')
end
if ~(isnumeric(Y))
   error('The second argument must be the output matrix.')
end
if ~(isreal(dt))
   error('The sampling time must be a positive real number.')
end
% if ~(isnegint(arg1) | (arg1==0) | ischar(arg1))
%    error('arg1 must be a non-positive integer or a string.')
% end
% if ~(isposint(arg2) | (arg2==0))
%    error('arg2 must be a non-negative integer.')
% end
if (num_args > 5)  % user specified approach
   if ~(strcmp(approach,'pseudo')|strcmp(approach,'corr')|strcmp(approach,'tvfil'))
      error('approach can only be ''pseudo'', ''corr'', or ''tvfil''.')
   end
else % user did not specify approach
   % Set default value for approach.
   approach = 'corr';
end
if (num_args > 6) % user specified confidence level
   if (strcmp(approach,'pseudo'))
      if isempty(conf_level)
         error('The confidence level must be a real number in the range (0 100) or must be set to NaN.')
      elseif ~((isreal(conf_level) & (conf_level>0) & (conf_level<100)) | isnan(conf_level))
         error('The confidence level must be a real number in the range (0 100) or must be set to NaN.')
      end      
   else % approach is not 'pseudo'
      if isempty(conf_level) 
         error('Confidence bounds are only computed when approach is ''pseudo''.')
      elseif ~(isnan(conf_level))
         error('Confidence bounds are only computed when approach is ''pseudo''.')
      end  
   end
else % user did not specify confidence level
   % Set default value for conf_level.
   conf_level = NaN;
end

% Set values for M1 and M2 based on arg1 and arg2.
if (isstr(arg1))
   hlength = arg2;
   if (strcmp(arg1, 'one'))
      M1 = 0;
      M2 = hlength-1;
   elseif (strcmp(arg1, 'two'))
      if (rem(hlength,2) == 0)
         M1 = -(hlength/2)+1;
         M2 = hlength/2;
      else
         M1 = - fix(hlength/2);
  	      M2 = fix(hlength/2);
      end
   else
      error('arg1 must be ''one'' or ''two'', or a non-positive integer.')
   end
else
      M1 = arg1;
      M2 = arg2;
      hlength = M2-M1+1; 
end

% Check that size of X is equal to size of Y.  If not, return error.
if ~((size(X,1)==size(Y,1)) & (size(X,2)==size(Y,2)))
   error('size(X) must equal size(Y)')
end

% Initialize some variables.
H_domain_start = M2+1;  
H_domain_end = size(X,1)+M1;
Hident = [];
bound = [];
sing_vectors = zeros(H_domain_end-H_domain_start+1,1);

% Remove the ensemble average from the input.
A = X' - ones(size(X,2),1)*mean(X');
X = A';
clear A

% Remove the ensemble average from the output. 
A = Y' - ones(size(Y,2),1)*mean(Y');
Y = A';
clear A

% Identify the time-varying system based on X and Y.
if (strcmp(approach,'tvfil'))  
   % use approach taken in tvfil
   
   cpu1 = cputime;


   for i = H_domain_start:H_domain_end

      % form the input matrix for time i
      Xi = X((i-M1):-1:(i-M2),:);
      Xi = Xi';
      % form the output vector
      yi = Y(i,:);
      yi = yi';
      % solve for hi using pinv function
      [pseudoinv,r] = pinv(Xi);
      hi = 1/dt*pseudoinv*yi;
      Hident = [Hident; hi'];

      % store the number of singular vectors used by pinv
      sing_vectors(i-H_domain_start+1) = r;
     
   end % for i
   
   cpu = cputime - cpu1;
   
elseif (strcmp(approach,'corr'))
   %find solution using correlation function approach but 
   % WITHOUT using Dave's technique to select singular vectors
   
   cpu1 = cputime;

   
   L = size(X,2);

   for i = H_domain_start:H_domain_end
      if (i == H_domain_start)
         % create the initial PHIxx matrix (use fact that it is symmetric)
         PHIxx = zeros((M2-M1+1),(M2-M1+1));
	      p = M2-M1+1;
	      for j = 1:(M2-M1+1)
	         q = p;
	         for k = j:(M2-M1+1)
	            PHIxx(j,k) = X(p,:)*X(q,:)';
	            PHIxx(k,j) = PHIxx(j,k);
	            q = q-1;
	         end
	         p = p-1;
	      end
	      PHIxx = 1/L*PHIxx;
      else
         % create PHIxx matrix using previous PHIxx matrix
	      newcov = [];
	      for k = M1:M2
	         newcov = [newcov X((i-M1),:)*X((i-k),:)'];
	      end
	      newcov = 1/L*newcov;
	 
	      newPHIxx = zeros((M2-M1+1),(M2-M1+1));
	      newPHIxx(2:(M2-M1+1),2:(M2-M1+1)) = PHIxx(1:(M2-M1),1:(M2-M1));
	      newPHIxx(1,:) = newcov;
	      newPHIxx(:,1) = newcov';
	      PHIxx = newPHIxx;		 
      end
   
      % create phiyx vector
      phiyx = [];
      for k = M1:M2
         phiyx = [phiyx; Y(i,:)*X((i-k),:)'];
      end
      phiyx = 1/L*phiyx;

      % find hi
      [pseudoinv,r] = pinv(PHIxx);
      hi = 1/dt*pseudoinv*phiyx;
      
      % append to Hident matrix 
      Hident = [Hident; hi'];

      % store the number of singular vectors used by pinv 
      sing_vectors(i-H_domain_start+1) = r;
                 
   end  % for i
   

   cpu = cputime - cpu1;
   
elseif (strcmp(approach,'pseudo'))
   %find solution using Dave's pseudoinv approach (adapted to tv case)
   
   cpu1 = cputime;

      
   L = size(X,2);
         
   for i = H_domain_start:H_domain_end
      if (i == H_domain_start)
         % create the initial PHIxx matrix (use fact that it is symmetric)
         PHIxx = zeros((M2-M1+1),(M2-M1+1));
	      p = M2-M1+1;
	      for j = 1:(M2-M1+1)
	         q = p;
	         for k = j:(M2-M1+1)
	            PHIxx(j,k) = X(p,:)*X(q,:)';
	            PHIxx(k,j) = PHIxx(j,k);
	            q = q-1;
	         end
	         p = p-1;
	      end
	      PHIxx = 1/L*PHIxx;
      else
         % create PHIxx matrix using previous PHIxx matrix
	      newcov = [];
	      for k = M1:M2
	         newcov = [newcov X((i-M1),:)*X((i-k),:)'];
	      end
	      newcov = 1/L*newcov;
	 
	      newPHIxx = zeros((M2-M1+1),(M2-M1+1));
	      newPHIxx(2:(M2-M1+1),2:(M2-M1+1)) = PHIxx(1:(M2-M1),1:(M2-M1));
	      newPHIxx(1,:) = newcov;
	      newPHIxx(:,1) = newcov';
	      PHIxx = newPHIxx;		 
      end
   
      % create phiyx vector
      phiyx = [];
      for k = M1:M2
         phiyx = [phiyx; Y(i,:)*X((i-k),:)'];
      end
      phiyx = 1/L*phiyx;

      % find the svd of PHIxx
      [U,S,V] = svd(PHIxx);
      
      % define tolerance to find rank of PHIxx (use default tolerance of rank() command, but don't use
      % rank() since it would call svd() once more and this would require a lot of computational 
      % effort for nothing)
      S = diag(S);
      tol = max(size(PHIxx)')*max(S)*eps;
      PHIxx_rank = sum(S>tol);
      
      % find Sinv
      Sinv = zeros(length(S),1);
      Sinv(1:PHIxx_rank) = 1./S(1:PHIxx_rank);
      Sinv = diag(Sinv);
      
      % compute output variance associated with each singular vector
      out_var = Sinv*((U'*phiyx).^2);
      
      % sort singular vectors and singular values according to 
      % corresponding output variance
      % [out_var, indices] = sort(out_var);
      % out_var = flipud(out_var);
      % indices = flipud(indices);
      % V = V(:,indices);
      % U = U(:,indices);
      % Sinv = diag(Sinv);
      % Sinv = Sinv(indices);
      % Sinv = diag(Sinv);
            
      % compute MDL for each candidate model
      cumul_out_var = cumsum(out_var);
      d = (1:hlength)';
      yvar = (std(Y(i,:)))^2;
      MDL = (ones(hlength,1)+log(L)/L*d).*(yvar*ones(hlength,1)-cumul_out_var);
      
      % compute pseudoinverse using optimal model order based on MDL
      [minMDL, r] = min(MDL);
      sing_vectors(i-H_domain_start+1) = r;
      pseudoinv = V(:,1:r)*Sinv(1:r,1:r)*U(:,1:r)';
      
      % solve for hi using the pseudoinverse
      hi = 1/dt*pseudoinv*phiyx;
      
      % append to Hident
      Hident = [Hident; hi'];
      
      if ~(isnan(conf_level))
         % compute the variance of the estimation error
         %resid_var = yvar - cumul_out_var(hlength);
         %estim_var = resid_var/L*diag(pseudoinv);
      
         % generate confidence bound for the estimation error
         % (assumes estimation error is normally distributed)
         %dom = [0:0.01:5];
         %gauss_pdf = exp(-(dom.^2)/2)/sqrt(2*pi);
         %gauss_cdf = cumsum(gauss_pdf)/sum(gauss_pdf);
         %boundary = dom(min(find(gauss_cdf>(conf_level/100))));
         %estim_bound = boundary*sqrt(estim_var);
      
         % generate confidence bound for total error
         % (for now, use the bound on the estimation error only)
         %boundi = estim_bound;       
         
         % append to bound matrix 
         %bound = [bound; boundi'];
      end % if
         
   end  % for i


   cpu = cputime - cpu1;
   
else
   
   error('No such approach.')
   
end  % if

% end tv_ident.m
