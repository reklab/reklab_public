function [coeffs,range,IRF,num_iterations] = tv_hammer_i(x,y,dt,order,M1,M2,tol)

% Identifies time-varying Hammerstein model between x and y using an
% iterative technique.  
%
% [coeffs,range,IRF,num_iterations] = tv_hammer_i(x,y,dt,order,M1,M2,tol)
%
% x:		input matrix
% y:		output matrix
% dt:		sampling interval
% order:	polynomial order
% M1:		- number of anticipatory lags in each IRF
% M2:		number of memory lags in each IRF
% tol:   stop iterating when VAFnew-VAFold < tol
%
% coeffs:	Cell array of polynomial coefficients describing the time-varying 
%		static nonlinearity.  Each row of the cell array contains a vector of
% 		polynomial coefficients corresponding to a specific instant in 
% 		time.
% range: 	Matrix containing the range of x used to identify each polynomial.  The 
%		ith row consists of the minimum and the maximum values of x used to 
%		identify the ith polynomial.
% IRF:		Matrix of weighting functions describing the time-varying linear dynamics.
%		The ith row consists of the weighting function for time i.
% num_iterations: Vector whose ith row indicates the number of iterations needed to 
%     identify the system at time i.
%

% Copyright 1999-2003, Robert E Kearney, Mireille Lortie
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


% Check that number of input arguments is valid.
num_args = nargin;
error(nargchk(7,7,num_args));

% Check that input arguments are valid.


% Check that size of x is equal to size of y.  If not, return error.
if ~((size(x,1)==size(y,1)) & (size(x,2)==size(y,2)))
   error('size(x) must equal size(y)')
end

% Initialize some variables.
[N,L] = size(x);
numlags = M2-M1+1;
IRF = zeros(N-numlags+1,numlags);
range = zeros(N-numlags+1,2);
for i = 1:(N-numlags+1)
   coeffs{i,1} = [];
end  % for i
num_iterations = zeros(N-numlags+1,1);

% Remove the ensemble average from the input.
A = x' - ones(size(x,2),1)*mean(x');
x = A';
clear A

% Estimate TV linear subsystem using output-input cross-correlation vector and 
% input autocorrelation matrix.

IRF = tv_ident(x,y,dt,M1,M2,'pseudo');% 27/10/08 TSV -changed from 'corr' to 'pseudo'

% For every time instant, use iterative technique to find the parameters of the static
% nonlinearity and of the linear subsystem that minimize the sum of squared differences
% between the observed and predicted outputs of the Hammerstein system at that time instant.
for i = (M2+1):(N+M1)
   
   %disp(['time instant: ' num2str(i-M2)])
   
   A = x((i-M1):-1:(i-M2),:);
   A = A';
   
   b = y(i,:);
   b = b';
   
   h = (IRF(i-M2,:))';

   VAFold = -Inf;
   keeplooping = 1;
   num_it = 0; 
   while (keeplooping)
            
      if ~(num_it == 0)
         VAFold = VAFnew;
         p_old = p;
      end
     
      if i == 501
          pause;
      end
      
      
      % Identify the static nonlinearity
      %disp('Identifying the static nonlinearity')
      C(:,1) = dt*sum(h)*ones(L,1);
      for j = 1:order
         C(:,j+1) = dt*(A.^j)*h;
      end % for j

%       TSV- replacing code below with built in regression function
      [pseudoinv,r] = pinv(C);
      p = pseudoinv*b;

%       [p,SE,pval,INmodel,stats,nextstep,history]=stepwisefit(C,b,'display','off');
      %remove parameters judged to not contribute much to fit
%       for k = 1:order
%           if INmodel(k) == 0
%               p(k) = 0;
%           end
%       end

      pred = C*p;
      p = p';
      p = fliplr(p);    %reformated for use with 'polyval'
      
      num_it = num_it + 0.5;

      % Compute VAF between observed and predicted outputs.  
      VAFnew = vaf(b,pred);
      %disp(['VAFnew = ' num2str(VAFnew)])
      
      if (VAFnew - VAFold) < tol
         if (VAFnew < VAFold)
            p = p_old;   
            num_it = num_it - 0.5;
         end
         break;
      end

      VAFold = VAFnew;
      h_old = h;
      
      % Identify linear subsystem
      %disp('Identifying the linear subsystem')
      mappedA = polyval(p,A);
      % DL and TSV - remove mean to remove offset from IRF - 0ct20/08
      mappedA = mappedA - (ones(length(mappedA),1)*mean(mappedA,1));
      [pseudoinv,r] = pinv(mappedA);
      h = (pseudoinv*b)/dt;
      pred = (mappedA*h)*dt;
      
      num_it = num_it + 0.5;
      
      % Compute VAF between observed and predicted outputs.  
      VAFnew = vaf(b,pred);
      %disp(['VAFnew = ' num2str(VAFnew)])
             
      if (VAFnew - VAFold) < tol
         if (VAFnew < VAFold)
            h = h_old;   
            num_it = num_it - 0.5;
         end
         break;
      end
   
   end  % while
   
   IRF(i-M2,:) = h';
   coeffs{i-M2,1} = p;
   range(i-M2,:) = [min(min(A)) max(max(A))];
   num_iterations(i-M2,1) = num_it;
   
end % for i


