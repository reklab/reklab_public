function vs = nl2vs (nl,vsin);
% Convert to  Volterra series from nl model  description
%
OrderMax = get(vsin,'vsOrderMax');

subsystems = get(nl,'elements');
p = subsystems{1};
p = nlident(p,'polyType','power');
mc=get(p,'polyCoef');

h = double(subsystems{2});
Ts = get(subsystems{2},'domainIncr');
hlen = length(h);

Q = get(p,'polyOrder');
if Q > OrderMax
    % truncate polynomial to aviod producing high-order kernels.
    str = ['truncating Volterra series at order ' num2str(OrderMax)];
    warning(str);
    Q = OrderMax;
end
kernels = cell(Q+1,1);


for q = 1:Q+1
    k = (mc(q)/Ts^(q-2))*hckern(h,q-1);
    kernels{q} = vkern(k,'kernOrder',q-1,'domainIncr',Ts);
end
vs=vsin;
set(vs,'elements',kernels,'vsOrderMax',Q,'nLags',hlen,...
    'comment','Transformed nlbl');
end
function kernel = hckern(g,Q);
% HCKERN - generates Q'th order Volterra kernel of a Hammerstein cascade.
% $Revision: 1.3 $
% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt

switch Q
    case 0
        kernel = sum(g);
    case 1
        kernel = g;
    case 2
        kernel = diag(g);
    otherwise
        % for kernels of order 3 or greater, we will need to use a multidimensinal
        % array.  First, create an empty kernel of the appropriate dimensions.
        
        hlen = length(g);
        indeces = hlen*ones(1,Q);
        kernel = zeros(indeces);
        
        % There should be a much more efficient way of doing this,
        command = 'kernel(i,i';
        idx = [',i'];
        for j = 3:Q
            command = [command, idx];
        end
        command = [command,') = g(i);'];
        for i = 1:hlen
            eval(command);
        end
end
end


