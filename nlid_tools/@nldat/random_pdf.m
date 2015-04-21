function y = random_pdf ( type, x, varargin) 
% usage: y = random_pdf ( type, x, varargin) 
% overlaid random_pdf for nldat objects
%  'normal',            mean     standard-deviation 
%  'chi-square',        dof
%  'student-t',         dof
%  'F-distribution'     dof1     dof2
%  'log-normal',        mean     variance   offset  
%  'rectangular',       xmin     xmax
%  'exponential',       mean
%  'binomial',          n        probabilty
%  'Poisson'            mean 

% Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

pdf_defs = { 'normal' 'chi_square' 'student-t' 'F-distribution' ...
      'log-normal' 'rectangular' 'exponential' 'binomial' ...
      'Poisson' };
options = { {'meanval' 0 'mean value' } ...
      {'sd' 1 'Standard Deviation' } ...
      {'dof' 20 'Degrees of Freedom' } ...
      {'dof2' 10 'Second degree of freedom' } ...
      {'offset' 0 'offset parameter' } ...
      {'minval' 0 'maximum value' } ...
      {'maxval' 1 'maximum value' } ...
      {'variance' 1 'variance'} ...
      {'n' 10 'number of samples'} ...
      {'prob' .1 'Probability of each event'} ... 
		};

arg_parse (options, varargin);
xd=double(x);
switch lower(type),
case 'normal'
   y=random_pdf(type,xd,meanval,sd);
case 'student-t'
   y=random_pdf(type,dof);
case 'F-distibution'
   y=random_pdf(type,dof, dof2);
case 'log_normal'
   y=random_pdf(type,meanval,variance,offset);
case 'rectangular'
   y=random_pdf(type,minval,maxval);
case 'exponential'
   y=random_pdf(type,meanval);
case 'biomial'
   y=random_pdf(type,n.prob);
case 'poisson'
   y=random_pdf(type,meanval);
otherwise
   disp(['Distribution not defied:' type]);
   disp(pdf_defs);
   error ('');
end
y=nldat(y);
set(y,'DomainValues',x,'Comment', [type ' pdf']);
% .... nldat/random_pdf






   

      
      
