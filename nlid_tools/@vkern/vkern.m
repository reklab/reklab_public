classdef vkern < kern
    % vkern - Volterra kernel class for NLID toolbox
    %  Parent: kern
    % Copyright 1999-2003, Robert E Kearney
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ..
    %
    
    properties
    end
    
    methods
        function VK = vkern (z, varargin)
            
            if nargin==0;
                return
            elseif nargin==1,
                VK=nlmkobj(VK,z); 
            elseif isa(z,'vkern')
                VK=nlmkobj(z,varargin{:})
            else
                VK=nlmkobj(VK,z,varargin{:});
            end
            return
        end

        
        function y = nlsim ( vk, u )
            % Simulate response of one Volterra kernel to an input signal
            % input options not fill defined as yet
            order =get (vk,'kernOrder');
            k=get(vk,'dataSet');
            Ts=get(vk,'domainIncr');
            ud=double(u);
            y=nldat(u);
            yd = vkfilt(k,ud,order)*Ts^order;
            set(y,'dataSet',yd);
            return
        end
    end
end
    function y = vkfilt(kern,u,order)
    % VKFILT - computes the output of an arbitrary-order Volterra kernel
    %
    %   y = vkfilt(kern,u);
    %
    %  kern is an "order" dimensional array of numbers containing the kernel.
    %  u is an N element vector containing the input.
    %  order is the order of the kernel.  If order is not specified, then the
    %  the number of matrix dimensions is used to determin the order.
    [urows,ucols] = size(u);
    u = u(:);
    y = zeros(size(u));
    N = length(u);
    
    if nargin < 3
        % figure out the order based on the number of kernel dimensions.
        nd = ndims(kern);
        if nd > 2
            order = nd;
        else
            % unfortunately, ndims returns 2 for scalars, vectors and matrices,
            % so we have to detect these separately.
            [nr,nc] = size(kern);
            if max(nr,nc)==1
                order = 0;
            elseif min(nr,nc)==1
                order = 1;
            else
                order = 2;
            end
        end
    end
    
    switch order
        case 0
            y = kern*ones(size(u));
        case 1
            y = filter(kern,1,u);
        case 2
            y = vk2filt(kern,u);
        otherwise
            ksize = size(kern);
            hlen = ksize(1);
            udel = u;
            for i = 1:hlen
                kk = kernel_slice(kern,i,order,hlen);
                y = y + vkfilt(kk,u,order-1).*udel;
                udel = [0;udel(1:N-1)];
            end
    end
    
    if ucols > 1
        y = y';
    end
    end
    
    
    
    
    function ks = kernel_slice(kern,slice,Q,hlen);
    % pull out a single Q-1 dimensional slice of an order Q kernel
    
    
    kcol = kern(:);
    numweights = hlen^(Q-1);
    start = (slice-1)*numweights+1;
    stop = slice*numweights;
    kscol = kcol(start:stop);
    command = 'ks = reshape(kscol,hlen';
    for i = 2:Q-1
        command = [command,',hlen'];
    end
    command = [command,');'];
    eval(command);
    end
    
    
    function y = vk2filt(k2,x);
    % compute the output from a second order kernel
    %
    %  syntax:  y=vk2filt(kern,x)
    %
    % dtw june 2003
    
    N = length(x);
    y = zeros(size(x));
    [nr,nc]=size(k2);
    
    % trikern is the triangular form of the kernel.
    % it is used to eliminate redundant computations,
    % since we can skip over mulitplications by zero.
    trikern = 2*triu(k2) - diag(diag(k2));
    
    xdel = x;
    for i = 1:nr
        y = y + filter(trikern(i,i:nc),1,xdel).*xdel;
        xdel = [0;xdel(1:N-1)];
    end
    end


