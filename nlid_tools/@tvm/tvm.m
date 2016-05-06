classdef tvm < nlm
    % define and create time-varying model objects
    %  Parent: nldat
    
    % Copyright 2000, Robert E Kearney
    % This file is part of the nlid toolbox, and is released under the GNU
    % General Public License For details, see ../copying.txt and ../gpl.txt
    properties
            tvFlag=true
    end
    
    methods
        
        function TVM=tvm(a,varargin);
            % Add object  specific parameters
            N.comment='IRF Model';
            if nargin==0;
                return
            elseif nargin==1,
                N=nlmkobj(N,a);
            elseif isa(a,'nlm')
                N=nlmkobj(a,varargin{:});
            else
                N.nlmkobj(N,a,varargin{:});
            end
        end
    end
end