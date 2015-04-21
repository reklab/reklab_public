function xAxisPanZomm(figNum)
% xAxisPanZoom - enable linked pan and zoom for linked axes 
%  
%   ... 
% Usage 
% xAxisPanZomm(figNum)
%
% figNum - optional argument specifying figure to use. Default is current
%           figure. 
% Example 
% zoo
% See also 
% 
%% AUTHOR: Robert Kearney 
%% $DATE: 04-Mar-2009 17:16:49 $ 
%% $Revision: 1.1 $ 
%% DEVELOPED: 7.7.0.471 (R2008b) 
%% FILENAME: xAxisPanZomm.m 
%Copyright 2002-2005 McGill University 
if nargin==1,
    figure(figNum);
end
fig = get(0,'CurrentFigure');
if isempty(fig), return; end
ax = findobj(fig,'Type','Axes');
nondatachild = logical([]);
for k=length(ax):-1:1
nondatachild(k) = isappdata(ax(k),'NonDataObject');
end
ax(nondatachild) = [];
linkaxes(ax,'x')
pan xon
zoom xon

