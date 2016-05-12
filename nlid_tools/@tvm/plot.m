function plot (M)
% plot function for tvm objects

% Copyright 2000, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../copying.txt and ../gpl.txt 

ModelType=get(M, 'Model_Type');
switch ModelType
case 'irf'
        Incr = get(M,'DomainIncr');
        Start=get(M,'DomainStart');
        dn=get(M,'DomainName');
        Z=get(M,'data');
        c=get(M,'comment');
       for i=1:length(Z),
           z(:,i)=get(Z{i},'data');
       end         
       % x axis is lag
       lag=domain(Z{i});
       lagname=get(Z{i},'domainname');
       % yaxis is time
       time=domain(M);
       m=mesh(lag,time,z');
       ylabel(dn);
       xlabel(lagname);
       
       title(c);
   case 'polynom'
       p=get(M,'data');
       nsamp=length(p);
       ds=get(M,'domainstart');
       di=get(M,'domainincr');
       t= ds + (0:nsamp-1)*di;
       range=get(p{1},'range');
       x=linspace(range(1),range(2));
       nreal=length(x);
       x1=repmat(x,nsamp,1);
       x1=reshape(x1,nsamp,1,nreal);
       x1=nldat(x1,'domainincr',di,'domainstart',ds);
       z=nlsim(M,x1);
       z=squeeze(double(z));
         m=mesh(x',t',z);
       xlabel('Input amplitude');
       ylabel('Domain');
       dn=get(M,'domainname'); ylabel(dn);
       c=get(M,'comment'); title(c);
   case 'nlbl'
       D=get (M,'data');
       for i=1:length(D),
           mi=D{i};
           sn{i,1}=mi{1,1};
           dl{i,1}=mi{1,2};
       end
       TVSN=M;
       set (TVSN,'model_type','polynom','data',sn);
       TVDL=M;
       set (TVDL,'model_type','irf','data',dl);
       subplot (1,2,1);
       plot (TVSN)
       subplot (1,2,2);
       plot (TVDL)
       c=get(M,'comment'); title(c);
 
      
   otherwise
       error (['tvm/plot - model type:'  ModelType ' not defined']);
   end
  
