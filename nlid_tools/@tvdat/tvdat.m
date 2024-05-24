classdef tvdat < nldat
    %tvdat Extension of nldat class to handle time-varying data
    %      sets such as spectorgrams
    
    properties
        realDomainIncr = 1 ;    % Realization domain increment
        realDomainName = 'Time {s}' ; % Realization domain name
        realDomainStart=0;   % Start value for realization domain
        realDomainValues=NaN;    % Real Domain Values - for use when realization are
                                 %sampled irregularly
    end
    
    methods
        function TV = tvdat (a, varargin)
            %tvdat  Construct an instance of this class
            if nargin==0,
                return
            elseif isa(a,'double');
                TV.dataSet=a;
                if nargin>1
                    set(TV,varargin{:});
                end
            elseif isa(a,'nldat');
                fields=fieldnames(a);
                nField=length(fields);
                for i=1:nField,
                    TV.(fields{i})=a.(fields{i});
                end
                if nargin > 1,
                    set (TV,varargin{:});
                end
            elseif isa(a,'tvdat');
                TV=nlmkobj(a,varargin{:});
            else
                set(TV, {a varargin{:}});
            end
        end
        
        function plot(TV)
            % Handle special case where these is only one realization
            [nSamp, nSig, nReal]=size(TV);
            if nReal==1,
                plot(nldat(TV));
                return
            end

            Z=squeeze(TV.dataSet);
            X=TV.domainValues;
            Y=TV.realDomainValues;
            waterfall(Y',X',Z);
            xlabel(TV.realDomainName);
            ylabel(TV.domainName); 
        end

         function mesh(TV)
            Z=squeeze(TV.dataSet);
            X=TV.domainValues;
            Y=TV.realDomainValues;
            mesh(Y',X',Z);
            xlabel(TV.realDomainName);
            ylabel(TV.domainName); 
            view(2);
            colorbar
        end
            
            
        
        
    end
end

