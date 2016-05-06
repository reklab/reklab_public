function yout = nlsim ( model, xin )
% Simulate response of a tvm 
% input options not fill defined as yet
% 

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

filter = get(model,'data');
if isa(xin,'double'),
    xin=nldat(xin);
    set(xin,'domainincr',get(model,'domainincr'));
end

x=double (xin);
incr = get (model,'domainincr');
P=get(model,'parameters');
assign(P);
%
% Simulate a time-varying polynomial response
%
ModelType=get(model,'Model_Type');
switch ModelType
case 'polynom'
    [nsamp, nchan, nreal]=size(xin);
    P=get(model,'data');
    for i= 1:nsamp;
        x=squeeze(xin(i,1,:));
        p=P{i};
        y=reshape(nlsim(p,x),1,1,nreal);
        if i==1,
            yout=y;
        else
            yout = cat (1,yout,y);
        end      
    end
    
    %
    % Simulate a time-varying IRF response
    %
case 'irf'
    delx = get(xin,'domainincr');
deli=get(model,'domainincr');
if delx ~= deli,
    W=(str2mat('Model & data have different domain increments', ...
        'the output of the IRF depends on the sampling rate', ...
        'Output may be scaled incorrectly and/or have the wrong increment'));
    warning(' ');disp(W)
end
    if NSides==1,
        sides='one';
    else
        sides='two';
    end
    x=x(:,1,:);  
    x=squeeze(x);
    for i=1:length(filter),
        f=double(filter{i});
        hirf(i,:)=f';
    end
    y = etvc(x,hirf,incr,sides);
    [n,m]=size(y);
    y=reshape(y,n,1,m);
    yout=xin;
    set(yout,'Comment','filtered','Data',y);
    %
    % Simulate a time-invariant response
    %
case 'nlbl'
    D=get(model,'data');
    for i=1:length(D),
        mi=D{i};
        p{i}=mi{1,1};
        l{i}=mi{1,2};
    end
    P=model;
    set(P,'Model_Type','polynom','data',p);
    L=model;
    set(L,'Model_Type','irf','data',l);
    y1=nlsim(P,xin);
    yout=nlsim(L,y1);
    
otherwise
    error ('Model_Type unknown');
end
