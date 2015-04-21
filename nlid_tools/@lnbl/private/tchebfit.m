function [coeff,vaf,out,mdls]=tchebfit(x,y, order,mode)
%
%	Tchebyschev polynomial fitting b/w input and output
%
%	Syntax: [coeff,vaf,out,mdls]=tchebfit(x,y,order, mode)
%
%		x	: input
%		y	: output
%		order	: polynomial order
%		mode 	: either 'auto', 'fixed' 'manual'
%				default is 'fixed'
%
%	in fixed mode, a polynomial of order order is estimated.  If the 
%	mode is either auto, or max, tchebyshev polynomials upto order
%	order are estimated.  In auto mode, a minimum description length
%	criterion are used to select the best model.  In manual mode, the 
%       model orders and vafs and mdl costs are displayed for the user to
%	use in selecting the best order.
%
%	The last three outputs are optional.  If specified, out will contain
%	the input, transformed by the computed Tchebyshev polynomial.

% Copyright 1991-2003, Robert E Kearney, David T Westwick and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 

x=x(:);
y=y(:);
b=max(x);
a=min(x);

%  rescale x to [-1,1]
x=(2*x - (b+a))/(b-a);
numcoeff = order + 1;
T=zeros(numcoeff,numcoeff);
A=zeros(length(x),numcoeff);
T(1,1)=1;
T(2,1:2)=[1 0];
A(:,1)=polyval(T(1,1),x);
A(:,2)=polyval(T(2,1:2),x);
for i=3:numcoeff
    T(i,1:i)=[2*T(i-1,1:(i-1)) 0] - [0 0 T(i-2,1:(i-2))];
    A(:,i)=polyval(T(i,1:i),x);
  end
if nargin == 3 
    mode = 'fixed';
  end
if mode(1) == 'f'
    coeff=A\y;
    out = A * coeff;
    coeff= [a;b;coeff];
    resid = y - out;
    vaf = 100 * (1 - ( std (resid) / std (y) )^2);
    mdls = mdl_cost(y,out,numcoeff,'y');
  else
    coeffs = [];
    vafs = [];
    mdls = [];

%      either manual or automatic polynomial order determination
%      calculate polynomial fits for all orders up to " order".
    for i = 0:order
	B = A(:,1:i+1);
	scratch = B\y;
	pad = zeros (order - i,1);
	scratch = [a;b;scratch;pad];
	coeffs = [coeffs, scratch];
	test = A * scratch(3:order+3);
	resid = y - test;
	vaf = 100 * (1 - ( std (resid) / std (y) )^2);
	vafs = [vafs, vaf];
        mdls = [mdls, mdl_cost(y,test,i,'y')]; 
      end
    if mode(1) == 'm'

%          Manual mode -- display variance accounted for by each
%          model order.  Let the user choose the model order.

	disp('   order      %VAF    MDL Cost')
	for i=0:order
	    disp ([i vafs(i+1) mdls(i+1)])
	  end
	ord = input ('Choose Polynomial Order: ');
	coeff = coeffs (1:ord+3,ord+1);
	vaf = vafs (ord+1);
      elseif mode(1) == 'a'

%          Automatic Mode.  Choose the order corresponding to the 
%          minimum in the MDL cost function

        [val,pos] = min(mdls);
        ord = pos - 1;
	coeff = coeffs (1:ord+3,ord+1);
	vaf = vafs (ord+1);
      end

    if nargout > 2
%
%             if requested by the user, calculate the polynomial output
%
	poly = [coeff(3:ord+3);zeros(numcoeff-ord-1,1)];
	out = A * poly;
      end    
  end

