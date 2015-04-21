function table_print (file,table,names,cw_in,form_in,prec_in)
% print table with heading
% Usage:
%
%   table_print (file,table,names,cw_in,form_in,prec_in)
%
%       file - output file (dump to screen if zero length)
%	table 
%	names 
%	cw_in integer vector - column width
%	form_in  - colum format ( integer(i), string(s),float(f))
%       prec_iin - column precision
%
%

% Copyright 2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 

if(length(names)>0),
	[m,n]=size(names);
else
	[n,m]=size(table);
end
if nargin <6,
	precision =[];
else
	precision = prec_in;
end
if  length(precision)==0,
	for i=1:m,
		precision(i)='3';
	end
end
if (nargin <5),
	form=[];
else
	form=form_in;
end
if length (form)==0,
	for i=1:m,
		form(i)='f';
	end
end
if (nargin <4),
	cw=[];
else
	cw=cw_in;
end
if length (cw)==0,
  [m,n]=size(names);
  n=max(n,8);
  for i=1:m,
    cw(i)=n+1;
  end
end
for i=1:80,
	blank_line(i)=' ';
end
header=blank_line;
current_col =0;
[m,n]=size(names);
for i=1:m,
   	name=deblank(names(i,:));
   	ilast=current_col+cw(i);
   	n1 = min(length(name),cw(i));
 	ifirst=ilast-n1+1;
	header(ifirst:ilast)=name(1:n1);
   	current_col=ilast+1;
    end
if(length(file)==0),
	disp(header);
else,
	fprintf(file,header);
	fprintf(file,'\n');
end
[nrow,ncol]=size(table);

for j=1:nrow,
   current_col=0;
   line=blank_line;
   for i=1:ncol,
	if( form(i)==abs('f')),
		fmt=['%' int2str(cw(i)) '.' precision(i) 'f'];
		s=sprintf(fmt,table(j,i));
      elseif(form(i)==abs('i')),
		s=int2str(table(j,i));
      else	
	s=sprintf(setstr(table(j,i)));
      end
      n1=length(s);
      if (n1 >cw(i)),
		s='*'
                n1=1;
	end
      ilast = current_col+cw(i);
      ifirst = ilast-n1+1;
      line(ifirst:ilast)=s(1:n1);
      current_col=ilast+1;
   end
if (length(file)==0),
   disp(line);
else,
	fprintf (file,line);
	fprintf (file,'\n');
end
end







