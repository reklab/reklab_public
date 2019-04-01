classdef pdf  < nldat
    % pdf - probabilty distribution function class for NLID toolbox
    
    properties
        parameterSet=param;
    end
    
    methods
        function sys = pdf  (a,varargin)
            % Add object  specific parameters
            sys.parameterSet(1)=param('paramName','pdfType',...
                'paramDefault','Density',...
                'paramHelp','Type of pdf',  ...
                'paramType','select',...
                'paramLimits',{'Density' 'Frequency' 'Probability','CDF' });
            sys.parameterSet(2)=param('paramName','nBins',...
                'paramDefault',100, ...
                'paramHelp','Number of bins; determined by p.domainValues if specified',  ...
                'paramType','number',...
                'paramLimits', [1 inf]);
            
            sys.comment='PDF ';
            if nargin==0;
                return
            elseif nargin==1,
                sys=nlmkobj(sys,a);
            elseif isa(a,'pdf')
                sys=nlmkobj(a,varargin{:});
            else
                sys=nlmkobj(sys,a,varargin{:});
            end
            
        end
        
        function p = nlident (p, z,  varargin)
            % CONSTRUCT an pdf  function object
            
            if nargin < 2,
                disp('NLIDtakes two inputs for pdf objects: pdf, z' );
            elseif nargin > 2,
                set(p,varargin);
            end
            assign (p.parameterSet);
            tempComment=p.comment;
            % PDF of data 
            if isa(z,'nldat') | isa (z,'double'),
                if isa(z,'double'),
                    z=nldat(z);
                end
                x=double(z);
                % Number of bins
                if ~isnan(p.domainValues(:))
                    set(p,'nBins',length(p.domainValues(:)));
                else
                    xMin=(min(x));
                    xMax=(max(x));
                    xRange=xMax-xMin;
                    if xRange==0,
                        xRange=1;
                    end
                        
                    delX=xRange/(nBins-1);
                    binCenters=(0:delX:xRange)+xMin;
                    p.domainValues=binCenters(:);
                end
                
                
                %
                pdfType=lower(pdfType);
                [dist,dom]=pdfnl(x,p.domainValues,pdfType);
                incr=dom(2)-dom(1);
                switch lower(pdfType)
                    case 'cdf'
                        set(p,'dataSet',dist(:),'chanNames',{'CumulativeProbabilty'}, ...
                            'chanUnits','Probabilty','comment','Cumulative Probabilty');
                    case 'density'
                        set(p,'dataSet',dist(:),'chanNames',{'Density'}, ...
                            'chanUnits','count','comment','Density');
                    case 'frequency'
                        set(p,'dataSet',dist','chanNames',{'Frequency'}, ...
                            'chanUnits','count','comment','Frequency');
                    case 'probability'
                        set(p,'dataSet',dist','chanNames',{'Probabilty'}, ...
                            'chanUnits',' ','comment','Probability');
                        
                end
                set (p,'comment',[ pdfType ' of ' char(z.chanNames)]);
                set (p,'domainName',[ 'Value of ' char(z.chanNames)]);
                set (p,'domainStart',min(dom), 'domainIncr',incr);;

                % pdf of randvar objects 
            elseif isa(z,'randvar'),
                pList= getParamValCell(z.parameterSet);
                if isnan(p.domainValues),
                    x=nlsim(z,ones(1000,1));
                    xMin=floor(min(x.dataSet));
                    xMax=ceil(max(x.dataSet));
                    xRange=xMax-xMin;
                    if xRange==0,
                        xRange=.1;
                    end
                    delX=xRange/(nBins-1);
                   xTemp=(0:delX:xRange)+xMin;
                   p.domainValues=xTemp(:);
     
                end
                     
                switch lower(pdfType)
                    case 'density'
                        switch lower(z.randvarType)
                            case 'uniform'
                                Y = unifpdf (p.domainValues, pList{:});
                            case 'chisquare'
                                Y = chi2pdf(p.domainValues, pList{:});
                           case 't'
                                Y = tpdf(p.domainValues, pList{:});
                            case 'f'
                                Y = fpdf(p.domainValues, pList{:});
                            case 'discrete uniform'
                                n=pList{1};
                                Y = unidpdf(p.domainValues, n);
                            otherwise
                                
              
                                PD = makedist(z.randvarType,pList{:});
                                
                                Y = pdf (PD, p.domainValues);
                        end
                        tempComment=['Probabilty density of randVar object'];
                    case 'cdf'
                        switch lower(z.randvarType)
                            case 'uniform'
                                Y = unifcdf (p.domainValues, pList{:});
                            case 'chisquare'
                                Y = chi2cdf(p.domainValues, pList{:});
                            case 'f'
                                Y = fcdf(p.domainValues, pList{:});
                            case 't'
                                Y = tcdf(p.domainValues, pList{:});
                            case 'discrete uniform'
                                Y = unidcdf(p.domainValues, pList{:});
                            otherwise
                                PD = makedist(z.randvarType,pList{:});
                                Y = cdf (PD, p.domainValues);
                        end
                        tempComment=['Cumulative Probabilty of randVar object'];
                    otherwise
                        error(['pdfType not supported for randVar objects:' pdfType]);
                end
                set(p,'dataSet',Y(:),'domainValues',p.domainValues(:));
            end
    
        end
        
        
        function plot ( C,varargin )
            % pdf/plotOverload plot function for PDF objects
            options={ {'help_flag' 0 'display help (0=No/1=yes)'} ...
                {'plotmode' 'bar' 'plotMode (bar/line/stem)'} ...
                {'linecolor' 'b' 'linecolor'} ...
                {'linewidth' 2 'Line width'} ...
                };
            if isstr(C),
                arg_help('nldat/plt',options);
                return
            end
            if arg_parse(options,varargin);
                return
            end
            
            assign(C.parameterSet);
            dx =C.domainIncr;
            c=C.dataSet;
            c = c(:);
            x = C.domainValues;
            nlc = length(linecolor);
            if nlc == 0
                linecolor = 'b';
            end
            if strcmp(plotmode, 'line')
                % C is an analytical PDF -- plot it normally
                plot(x,c,linecolor);
                ylabel(ptfType);
            elseif strcmp(plotmode, 'stem')
                % C is an analytical PDF -- plot it normally
                h=stem(x,c);
                set(h,'color',linecolor, 'linewidth',linewidth)
                ylabel(pdfType);
            elseif strcmp(plotmode,'add')
                h=line(x,c);
                set(h,'color',linecolor,'linewidth',linewidth);
            else
                % C is an experimental PDF -- plot it like a histogram
                x2 = [x];
                c2 = [c];
                bar(x2,c2,1,linecolor);
                set (gca,'xlim',[ min(x)-dx max(x)+dx],'ylim',[ 0 (max(c2))]);
                ylabel(pdfType);
            end
            
            xlabel('Value')
            title(C.comment);
        end
        
        
    end
    
    
    
end
function [dist,dom] = pdfnl(x,binCenters,option)
% compute the probabilty distibution of a signal
%	THE FUNCTION PDF COMPUTES THE PROBABILITY DISTRIBUTION FUNCTION
%	OF A GIVE CHANNEL X.
%
%	USAGE:	[dist,dom] = pdf(x,numbins,option)
%
%	x	: input vector
%	numbins: number of bins used in order to calculate the
%		  discrete PDF or normalized histogram.
%	option  - string determined result
%                     'freq' - frequency histogram Ni
%                     'prob' - probabilty function Ni/Ntotal
%                     'dens' - proabilty density Ni/Ntotal*BinWidth
%                     'cdf - cumulativeProbabilty'


% Copyright 1991-2003, Robert E Kearney and Eric J Perreault
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see ../copying.txt and ../gpl.txt

if (nargin==2),
    option='probability';
end
%% Get rid of NANs
iNan=find(~isnan(x));
xtemp=x(iNan);
x=xtemp;
%
option=lower(option);
scale = length(x);
[dist,dom]=hist(x,binCenters);
dist=hist(x,dom);
if (strcmp(option,'frequency')),
    dist = dist;
elseif(strcmp(option,'cdf')),
    dist =cumsum( dist / scale);
elseif(strcmp(option,'probability')),
    dist = dist / scale;
elseif(strcmp(option,'density')),
    binwidth=unique(chop(diff(dom),5));
    if length(binwidth)>1,
b        error('Different bindwidths not yet supported');
    end
    if dist(1) ~=0 | dist(end)~=0,
        % warning('first and last bins must be empty for density estimates');
    end
    dist = dist / (binwidth*scale);
else
    error(['pdf - bad option:' option]);
    
end
end



