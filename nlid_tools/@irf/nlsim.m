function y = nlsim ( model, xin )
% irf/nlsim Simulate response of IRF to input data set
% input options not fill defined as yet
filter = model.dataSet;
if isa(xin,'double'),
    xin=nldat(xin);
    set(xin,'domainIncr',model.domainIncr);
end
delx = xin.domainIncr;
deli=model.domainIncr;
if delx ~= deli,
    W=(str2mat('Model & data have different domain increments', ...
        'the output of the IRF depends on the sampling rate', ...
        'Output may be scaled incorrectly and/or have the wrong increment'));
    warning(' ');disp(W)
end
x=double (xin);

incr = model.domainIncr;
assign(model.parameterSet);
%
% Simulate a time-varying response
%
[irfLen, irfDim, nSampIrf]=size(model);
 if nSides==1,
      offSet=0
      offSetStart=irfLen-1;
   else
      offsetStart=(irfLen-1)/2;
      offsetEnd=-offSetStart;
   end
[nSamp, nChan, nReal]=size(xin)
if (tvFlag),
  
   for iReal=1:nReal,
       for iSamp=1:nSamp,
           curIRF=model(:,1,iSamp);
           iEnd=iSamp+offSetEnd;
           iStart=iSamp-offSetStart;
       end
   end
   x=x(:,1,:);  
   x=squeeze(x);
   filter=squeeze(filter);
   filter=filter;
   yout = etvc(x,filter',incr,sides);
   [n,m]=size(yout);
   yout=reshape(yout,n,1,m);
   y=xin;
   set(y,'c','filtered','Data',yout);
   %
   % Simulate a time-invariant response
   %
else
  x=x(:,1,:);  
    
    [nsamp, nchan,nreal]= size(filter);
    for i=1:nchan,
        for j=1:nreal,
            yout(:,i,j) = filter_ts(filter(:,i,j), x(:,i,j), nSides, incr);
        end
    end
    y=xin;
    set(y,'comment','filtered','dataSet',yout);
    
end
set (y,'chanNames',{'Predicted output'});


