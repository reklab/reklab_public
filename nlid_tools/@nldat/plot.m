function plot (d, varargin)
% nldat/plt - overloaded plot function for "nldat" class
% plot (d, varagin);
% d - data to plot
% varargin - option/value pairs use  plot(d,'?' ) -to show options
% options={{'plotmode' 'line' 'plot type (line/xy/Super/bar/stem)'} ...
%         {'help_flag' 0 'display help (0=No/1=yes)'} ...
%           {'line_color' '' 'Line color'} ...
%            {'LineWidth' 0.5 'Line color'} ...
%         {'nh' NaN 'Number of horizontal plots'} ...
%         {'nv' NaN 'Number of vertical plots'} ...
%         {'nplt' 1 'Plot Number' } ...
%         {'realizationMode' 'super' 'mode for plotting realizations (super/offset/mesh/waterfall)' } ...
%          {'realizationOffset' 1 'offset for plotting realizations' } ...
%         { 'xpanzoom' false 'enable linked zoom and pan for all axes'} ...
%         { 'xpanwidth' NaN 'width of pan window' } ...
%         {'xmode' 'linear' 'x axes mode [linear/log' } ...
%         {'ymode' 'linear' 'y axes mode [linear/log/db' } ...
%
%
% % Copyright 1999-2003, Robert E Kearney
% This file is part of the nlid toolbox, and is released under the GNU
% General Public License For details, see copying.txt and gpl.txt

options={{'plotmode' 'line' 'plot type (line/xy/Super/bar/stem)'} ...
    {'help_flag' 0 'display help (0=No/1=yes)'} ...
    {'lineColor' '' 'Line color'} ...
    {'LineWidth' 0.5 'Line color'} ...
    {'nh' NaN 'Number of horizontal plots'} ...
    {'nv' NaN 'Number of vertical plots'} ...
    {'nplt' 1 'Plot Number' } ...
    {'realizationMode' 'super' 'mode for plotting realizations (super/offset/mesh/waterfall)' } ...
    {'realizationOffset' 1 'offset for plotting realizations' } ...
    { 'xpanzoom' false 'enable linked zoom and pan for all axes'} ...
    { 'xpanwidth' NaN 'width of pan window' } ...
    {'xmode' 'linear' 'x axes mode [linear/log' } ...
    {'ymode' 'linear' 'y axes mode [linear/log/db' } ...
    
    };
if nargin==2 & strcmp(varargin{1},'?'),
    arg_help('nldat/plt',options);
    return
end
if arg_parse_c('exact',options,varargin);
    return
end
%
% Determine default layout
%
[nsamp, nchan, nreal]= size(d);
if isnan (nv) & isnan (nh),
    % neither nv or nh set
    nh=ceil(nchan/2);
    nv=ceil(nchan/nh);
elseif isnan (nv),
    % nh set
    nv=ceil(nchan/nh);
elseif isnan(nh),
    % nv set
    nh = ceil (nchan/nv);
end
;

incr=d.domainIncr;
start=d.domainStart;
names=d.chanNames;
if isnan(d.domainValues),
    t= ((1:nsamp)-1)*incr +start;
else
    t=d.domainValues;
end
if strcmp(xmode,'log'),
    t=log10(t);
end
x=double(d);
if strcmp(ymode,'log'),
    x=log10(x);
end
if strcmp(ymode,'db'),
    x=20*log10(x);
end

%
% xy plot
%
plotmode=lower(plotmode);
if strcmp(plotmode,'xy'),
    for i=1:2:nchan,
        plot(x(:,i,:),x(:,i+1,:),'.')
    end
    xlabel(d.chanNames{1});
    ylabel(d.chanNames{2});
        title(d.comment);
elseif strcmp(plotmode,'super'),
    plot (t,x);
    xlabel(d.domainName);
    yunits = d.chanUnits;
    if ~isempty(yunits)
        ylabel(d.chanUnits);
    end
    title (d.comment);
    
    % this next hack produces a legend with all of the channel names
    % harvested from the ChanNames cell array.
    numchan = size(x,2);
    cmd = 'legend(''';
    for i = 1:numchan
        cmd = [cmd  d.chanNames{i} ''','''];
    end
    cmd = [cmd(1:end-1) '''NorthEast'')'];
    eval(cmd);
    % end of legend hack.
        title(d.comment);
    
    
elseif strcmp(plotmode,'bar'),
    bar (t,x,1);
        title(d.comment);
elseif strcmp(plotmode,'stem'),
    stem(t,x);
        title(d.comment);
else
    %
    % Line plots
    %
    if isnan(nplt),
        nplt=1;
    end
    
    for i=1:nchan,
        if (nchan > 1)
            NP=nplt+i-1;
            subplot (nv,nh,NP);
        end
        xd = squeeze(x(:,i,:));
        if nreal==1,
            nlc = length(lineColor);
            if nlc>0,
                plot (t,xd,lineColor(max(nlc,i)));
            else 
                plot (t,xd,'linewidth',LineWidth);
            end
            
            if strcmp(xmode,'log'),
                xlabel([ 'Log ' d.domainName]);
            else
                xlabel(d.domainName);
            end
            switch ymode
                case 'log'
                    ylabel([' log (' names{i} ')']);
                case 'db'
                    ylabel([names{i} ' dB']);
                otherwise
                    ylabel(names{i});
            end
            
            
            
        elseif nreal > 1,
            y=(1:nreal);
            if strcmp(realizationMode,'surf'),
                surf(y,t,xd);
            elseif strcmp(realizationMode,'mesh');
                mesh(y,t,xd);
            elseif strcmp(realizationMode,'waterfall');
                waterfall(y,t,xd)
            elseif strcmp(realizationMode,'offset');
                for ir=1:nreal;
                    xd(:,ir)=xd(:,ir)+realizationOffset*(ir-1);
                end
                
                plot (t,xd);
            elseif strcmp(realizationMode,'super');
                nlc = length(lineColor);
                if nlc>0,
                    plot (t,xd,lineColor(max(nlc,i)));
                else
                    plot (t,xd);
                end
            end
        end
        if i==1,
            title(d.comment);
        end
        
        
        
    end
    if xpanzoom,
        xAxisPanZoom (gcf);
    end
    if ~isnan(xpanwidth),
        xStart=get(d,'domainStart');
        xEnd=xStart+xpanwidth;
        set (gca,'xlim', [ xStart xEnd]);
    end

    
end


end
%	.../@nldat/plot
