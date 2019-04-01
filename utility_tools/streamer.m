function handle = streamer(TitleString,yPos)
%  STREAMER  Titles for an entire figure.
% 	     STREAMER('text', yPos) adds text at the top of the current figure,
% 	     going across subplots. 
% [ypos] - optional postion of text
% 
% 	     See also XLABEL, YLABEL, ZLABEL, TEXT, TITLE.
%
%	Keith Rogers 11/30/93
% 22 Sep 94 REK Modify for 4.2
% 26 Aug 97 REK Modify for 5.1
% 26 Aug 98 REK Problem with initializing StreamerHand
if nargin<2,
    yPos=0.92;
end
StreamerHand=[];
ax = gca;
fig = gcf;
sibs = get(fig, 'Children');
for i = 1:max(size(sibs))
	if(strcmp(get(sibs(i),'Type'),'axes'))
		sibpos = get(sibs(i),'Position');
		if (sibpos(2)+sibpos(4) == .93)
				StreamerHand = sibs(i);
		end
	end
end
if (isempty(StreamerHand)),
	StreamerHand = axes('position',[.1 yPos .8 .05],'Box','off','Visible','off');
else
	axes(StreamerHand);
end
title(TitleString,'visible','on');
handle = get(gca,'title');
axes(ax);



