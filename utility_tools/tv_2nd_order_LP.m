
function [H] = tv_2nd_order_LP(I,B,K,hlength,dt,dt1)

% Generates matrix of second-order low-pass filters (IRF's). 
%
% Usage: [H] = tv_2nd_order_LP(I,B,K,hlength,dt,dt1)
%
% I,B,K:      vectors containing values of the inertial, viscous, and elastic parameters
%             (each vector element corresponds to a specific instant in time)
% hlength:    length of the IRF at every point in time
% dt:         sampling interval
% dt1:        sampling interval initially used to sample the continuous-time 2nd-order IRF
%             If dt1 is not the equal to dt, the continuous-time IRF is first sampled at a 
%             sampling freq. of 1/dt1, then it is anti-alias filtered and resampled at a 
%             sampling freq. of 1/dt.  This avoids aliasing effects associated with sampling the
%             continuous-time IRF.  (By default, dt1=dt.  In such a case, the user must make sure
%             that 1/dt is MUCH larger than the desired cutoff, so that aliasing effects are minimal.)
%             If dt1 is specified, dt must be an integer multiple of dt1.
%
% H:          matrix of second-order low-pass filters (IRF's)
%
% The matrix H has the following format:
%     H:  ==== Lag Time ====>
%               ||
%               ||
%               ||
%            Discrete
%              Time
%               ||
%               ||
%               ||
%               ||
%               \/
%

% Check that number of input arguments is valid.
num_args = nargin;


% Check that input arguments are valid.
if ~(isnumeric(I) & (ndims(I)==2) & ((size(I,1)==1) | (size(I,2)==1)))
   error('I must be a column or a row vector.')
end
if ~(isnumeric(B) & (ndims(B)==2) & ((size(B,1)==1) | (size(B,2)==1)))
   error('B must be a column or a row vector.')
end
if ~(isnumeric(K) & (ndims(K)==2) & ((size(K,1)==1) | (size(K,2)==1)))
   error('K must be a column or a row vector.')
end
if ~(isposint(hlength))
   error('hlength must be a positive integer.')
end
if ~(isposreal(dt))
   error('dt must be a positive real number.')
end
if (num_args==6)
   if ~(isposreal(dt1))
      error('The sampling interval used to sample the continuous-time IRF must be a positive real number.')
   end
else
   dt1 = dt;
end

% Make sure that I, B, and K are of the same length.
if ~((length(B)==length(I)) & (length(K)==length(I)))
   error('I, B, an K must have the same length.')
end

% Make sure that dt is an integer multiple of dt1.
if ~(dt1==dt)
   ratio = dt/dt1;
   if ~(isposint(ratio))
      error('The sampling interval dt must be a multiple of the sampling interval used to sample the continuous-time IRF.');
   end
else
   ratio = 1;
end

% Compute lag axis.
hlength1 = hlength*ratio;
lag = 0:dt1:(hlength1-1)*dt1;

% Generate the matrix of IRF's
H = [];
for i = 1:length(I)
   h = fkbi(lag, [K(i) B(i) I(i)]);
   if ~(dt1==dt)
      h = resample(h,1,ratio);
   end
   H = [H; h'];
end

% end tv_2nd_order_LP.m








