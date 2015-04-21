classdef kern < nldat
    
    % kern - kernel class for NLID toolbox 
    %  Parent: nldat
    
    properties
       parameterSet
     
    end

    
    methods
        function K = kern (x,varargin)
     
            K.parameterSet=param('paramName','nLags','paramDefault',16, ...
                'paramHelp','Number of lags', ...
                'paramType','number','paramLimits', { 1 inf});
            K.parameterSet(2)=param('paramName','nSides','paramDefault',1,...
                'paramHelp','number of sides', ...
                'paramType','number','paramLimits',{1 2});
            K.parameterSet(3)=param('paramName','kernOrder','paramDefault',1,...
                'paramHelp','kernel order', ...
                'paramType','number','paramLimits',{1 2});
            K.parameterSet(4)=param('paramName','tvFlag','paramDefault','No', ...
                'paramHelp','tvFlag(Yes/No)');
         
            if nargin ==0,
                return
            elseif isa(x,'kern'),
                K=nlmkobj(x,varargin{:});
            else
                K=nlmkobj(K,x,varargin{:});
            end
        end
        
        function hh = plot ( C )
            % Overloaded plot function for kernel objects
            assign(C.parameterSet);
           
            Ts = C.domainIncr;
            c=C.dataSet;
            t=domain(C);
            tlim=[ min(t) max(t)];
            if tlim(1)==tlim(2),
                tlim(2)=tlim(1)+1;
            end
            switch kernOrder
                case 0
                    handle = stem(c,'filled');
                    set(gca,'xticklabel',[]);
                    ylabel('zero-order kernel');
                case 1
                    handle=gca;
                    plot(nldat(C));
%                     handle = plot (t,c);
%                     set(gca,'xlim',tlim);
%                     ylabel(C.chanNames);
%                     xlabel(C.domainName);
%                     title(C.comment);
                case 2
                    handle = surface (t,t,c); view (45,30);
                    set (gca,'xlim',tlim,'ylim',tlim);
                    colormap('jet');caxis('auto');
                    shading ('interp');
                    zlabel('second-order kernel');
                    xlabel('lag (sec.)');
                    ylabel('lag (sec.)');
                case 3
                    [xx,yy,zz] = meshgrid(t,t,t);
                    hlen = length(t);
                    zslices = t([1,floor(hlen/3),floor(2*hlen/3),hlen]);
                    xslice = t(ceil(hlen/2));
                    handle = slice(xx,yy,zz,c,[],[],zslices);
                    set(handle,'edgecolor','none','facecolor','interp',...
                        'facealpha','interp');
                    alpha('color');
                    alphamap('default');
                    alphamap('increase',0.1);
                    caxis('auto');
                    colormap('jet');
                    axis tight
                    xlabel('lag (sec.)');
                    ylabel('lag (sec.)');
                    zlabel('lag (sec.)');
                    cb = colorbar('vert');
                    axes(cb)
                    ylabel('third-order kernel');
            end
            if nargout == 1
                hh = handle;
            end
        end
        
        function k = nlident ( k, z, varargin )
            % Default identification for kernel objects
            % simply set the kernel values equal to z
            z=nldat(z);
            k = mvField(fieldnames(nldat), z,k);
            if length(varargin) > 0,
                set (k,varargin{:});
            end
            
        end
        
        function ks = smo(k,n)
            % overloaded smooth function for kern objects
            % syntax ks = smo(k,n);
            % n is the number of passes (default 1).
            if nargin<2,
                n=1;
            end
            ks=k;
            kdat=k.dataSet;
            order =get(k,'kernOrder');
            switch order
                case 0
                    % zero order kenel -- do nothing
                    kdat = kdat;
                case 1
                    % first order kernel -- smooth it.
                    kdat = smo(kdat,n);
                case 2
                    % second order kernel -- use a 2d smoother
                    kdat = smo_2d(kdat,n);
                case 3
                    [nr,nc,nl] = size(kdat);
                    for i = 1:nl
                        kdat(:,:,i) = smo_2d(squeeze(kdat(:,:,i)),n);
                    end
                    for i = 1:nc
                        kdat(:,i,:) = smo_2d(squeeze(kdat(:,i,:)),n);
                    end
                    for i = 1:nr
                        kdat(i,:,:) = smo_2d(squeeze(kdat(i,:,:)),n);
                    end
                otherwise
                    error('smo not implemented for kernels of order 4 and up');
            end
            cc = k.comment;
            if strcmp(cc,'Default comment')
                cc = 'smoothed';
            else
                cc = ['smoothed ' cc];
            end
            set(ks,'dataSet',kdat,'comment',cc);
        end
        
        
        
    end
end
function C = mvField (fieldList, A,B)
% Move  fields from A to B
C=B;
for i=1:length(fieldList),
    curName=fieldList{i};
    C.(curName)= A.(curName);
end
end

function ks = smo_2d(k,n)
% smo will smooth a matrix column by column -- so we have to transpose
% to smooth all of the rows, and transpose back again.
ks = smo(k,n);
ks = smo(ks',n)';
end




