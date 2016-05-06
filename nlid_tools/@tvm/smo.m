function y=smo(x,N);
% smo for tvm objects

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

MT=get(x,'Model_Type');
D=get(x,'data');
L=length(D);
switch MT
case 'polynom'
case 'irf'
    for i=1:L
        dtemp=D{i};
        dtemp=smo(dtemp,N);
        D{i}=dtemp;
    end
    y=x;
    set(y,'data',D);
case 'nlbl'
    for i=1:L
        dtemp=D{i};
        el=get(dtemp,'elements');
        Itemp=el{1,2};
        Itemp=smo(Itemp,N);
        el{1,2}=Itemp;
        set(dtemp,'elements',el);
        D{i}=dtemp;
    end
    y=x;
    set(y,'data',D);
    
otherwise
    error ('Model Type Note supported');
end
