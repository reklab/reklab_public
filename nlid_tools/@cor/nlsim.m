function y = nlsim ( model, xin )
% Simulate response of IRF to input data set
% input options not fully defined as yet

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

filter = model{1};
xin=nldat(xin);
x=double (xin);
x=x(:,1);
incr = get (model,'Ts');
numsides=get(model,'NSides');
[nsamp, nchan,nreal]= size(filter);
for i=1:nchan,
   for j=1:nreal,
      yout(:,i,j) = filter_ts(filter(:,i,j), x, numsides, incr);
   end
end
y=xin;
set(y,'Comment','filtered','Data',yout);
