function mm(action)
%MM  Memory Monitor displays runtime memory information
%   MM(t) updates data in t seconds steps
%
%   version v1.4
%   author  Elmar Tarajan [MCommander@gmx.de]
%
if eval(version('-release'))<14 | ~ispc
   error('Sorry...MATLAB R14 on Windows PC is required. :(')
end% if
%
if nargin==0 | isa(action,'double')
   %
   if ~isempty(findall(0,'Name','mm v1.4'))
      disp('warning: only one instance is allowed');
      return
   end% if
   check = 'off';
   if ~exist('action','var')
      action = 5;
      check = 'on';
   end% if
   %
   ftr.units = 'pixels';
   figure(ftr,'menubar','none','Name','mm v1.4', ...
      'NumberTitle','off','CloseRequestFcn','mm(''X'')', ...
      'position',[4 get(0,'ScreenSize')*[0;0;0;1]-80 137 80], ...
      'resize','off','UserData',action)
   ftr.style = 'edit';
   ftr.enable = 'off';
   ftr.backgroundcolor = [.635 .635 .635];
   uicontrol(ftr,'position',[1 19 137 61])
   uicontrol(ftr,'position',[1  1  69 18])
   uicontrol(ftr,'position',[70 1  68 18])
   %
   ftr.style = 'text';
   ftr.enable = 'on';
   ftr.backgroundcolor = [.784 .78 .6];
   ftr.foregroundcolor = [.294 .294 .294];
   uicontrol(ftr,'position',[4 64 65 13],'String','RAM')
   uicontrol(ftr,'position',[4 50 65 13],'String','MATLAB')
   uicontrol(ftr,'position',[4 36 65 13],'String','JAVA')
   uicontrol(ftr,'position',[4 22 65 13],'String','SWAP')
   uicontrol(ftr,'position',[3  4 31 12],'Tag','inmem', ...
      'tooltipstring','Number of M-Files in memory (inmem)')
   uicontrol(ftr,'position',[36 4 31 12],'Tag','hg', ...
      'tooltipstring','Number of handle graphics objects in memory')
   %
   ftr.backgroundcolor = [.804 .839 .894];
   ftr.fontname = 'small fonts';
   ftr.fontsize = 5;
   ftr.HorizontalAlignment = 'right';
   uicontrol(ftr,'position',[70 22 65 55],'string','V1.4 ')
   %
   ftr.fontname = 'small fonts';
   ftr.fontsize = 7;
   ftr.HorizontalAlignment = 'left';
   ftr.foregroundcolor = [.9 .9 .9];
   ftr.backgroundcolor = [.392 .561 .702];
   uicontrol(ftr,'position',[71 65 50 11],'Tag','RAM')
   uicontrol(ftr,'position',[71 51 50 11],'Tag','MATLAB')
   uicontrol(ftr,'position',[71 37 50 11],'Tag','JAVA')
   uicontrol(ftr,'position',[71 23 50 11],'Tag','SWAP')
   %
   ftr.backgroundcolor = [.659 .753 .835];
   uicontrol(ftr,'position',[114 51 20 11],'Tag','LARGE')
   %
   ftr.callback = 'mm(get(gcbo,''String''))';
   ftr.fontname = 'MS Sans Serif';
   ftr.fontweight = 'bold';
   ftr.backgroundcolor = [.8 .8 .8];
   ftr.foregroundcolor = [.2 .2 .2];%[1 .573 .141];
   ftr.style = 'pushbutton';
   uicontrol(ftr,'position',[121 3 15 14],'string','X', ...
      'tooltipstring',['Matlab Memory Monitor v1.4' char(10) ...
      'by Elmar Tarajan [MCommander@gmx.de]'])
   ftr.style = 'pushbutton';
   %   ftr.foregroundcolor = [1 .573 .141];
   ftr.uicontextmenu = uicontextmenu('parent',gcf);
   uicontrol(ftr,'position',[73 3 48 14],'string','pause','UserData',action)
   uimenu(ftr.uicontextmenu,'Label', sprintf('Timer (%.1fsec)',action))
   uimenu(ftr.uicontextmenu,'Label','1','separator','on','Callback','mm(''periode'')')
   uimenu(ftr.uicontextmenu,'Label','2','Callback','mm(''periode'')')
   uimenu(ftr.uicontextmenu,'Label','5','Callback','mm(''periode'')','checked',check)
   uimenu(ftr.uicontextmenu,'Label','10','Callback','mm(''periode'')')
   uimenu(ftr.uicontextmenu,'Label','20','Callback','mm(''periode'')')
   set(gcf,'Handlevisibility','off')
   drawnow
   mm('start')
   return
end% if
%
try
   tmp = set(0,'HandleVisibility');
   set(0,'HandleVisibility','on');
   %
   switch action
      case 'pause'
         T = timerfind('Name','mm');
         stop(T)
         delete(T)
         for i={'RAM' 'MATLAB' 'JAVA' 'SWAP' 'inmem' 'hg'}
            set(findall(0,'Tag',char(i)),'Enable','off')
         end% for
         set([gcbo findall(0,'String','X')],'Enable','on')
         set(gcbo,'String','start')% ,'foregroundcolor',[ 0.2 0.5 0])
         %
      case 'start'
         set(gcbo,'String','pause')% ,'foregroundcolor',[1 .573 .141])
         for i={'RAM' 'MATLAB' 'JAVA' 'SWAP' 'inmem' 'hg'}
            set(findall(0,'Tag',char(i)),'Enable','on')
         end% for
         T = timer('Name','mm', ...
            'TimerFcn','mm(''memory'')', ...
            'Period',get(findall(0,'Name','mm v1.4'),'UserData'), ...
            'ExecutionMode','fixedDelay');
         start(T);
         %
      case 'periode'
         tmp = findobj(gcbf,'type','uimenu');
         t = get(gcbo,'Label');
         set(tmp(end),'Label',['Timer (' t 'sec)'])
         set(findobj(gcbf,'type','uimenu'),'checked','off')
         set(gcbo,'checked','on')
         T = timerfind('Name','mm');
         if ~isempty(T)
            stop(T)
            set(T,'Period',str2num(t))
            start(T)
         else
            set(findobj(0,'Name','mm v1.4'),'UserData',str2num(t))
         end% if
         %
      case 'X'
         try
            mm('pause')
         end% try
         closereq
         %
      case 'memory'
         %
         % inmem
         set(findall(0,'Tag','inmem'),'string',length(inmem))
         %
         % hg
         set(findall(0,'Tag','hg'),'string',length(findall(0)))
         %
         % JAVA
         total = java.lang.Runtime.getRuntime.totalMemory;
         inuse = total-java.lang.Runtime.getRuntime.freeMemory;
         update('JAVA',inuse/1048576,total/1048576)
         %
         % MATLAB
         mem = dataread('string',evalc('feature(''memstats'')'),'%s');
         mem = mem(find(strcmp(mem,'MB'))-1);
         mem = eval(['[' sprintf('%s ',mem{:}) ']']);
         %
         tmp = evalin('base','whos');
         inuse = sum([tmp.bytes])/1048576;
         total = mem(end)+inuse;
         %
         update('MATLAB',inuse,total)
         h = findall(0,'Tag','LARGE');
         pos = get(findall(0,'Tag','MATLAB'),'position');
         set(h,'ToolTipString',sprintf('Largest available memory block is %dMB',mem(10)), ...
            'position',get(h,'position').*[0 1 0 1]+[pos(1)+pos(3) 0 mem(10)/total*49 0])
         %
         % RAM
         update('RAM',mem(1),mem(3))
         %
         % SWAP
         update('SWAP',mem(4),mem(6))
         %
   end% switch
   set(0,'HandleVisibility','off');
end% try
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update(tag,inuse,total)
%-------------------------------------------------------------------------------
tmp = inuse/total;
h = findall(0,'Tag',tag);
info = sprintf('in use->%.1fMB|%.1fMB<-total',inuse,total);
set(findall(0,'String',tag),'tooltipstring',info);
set(h,'string',sprintf('%.0f%%',tmp*100),'tooltipstring',info, ...
   'position',get(h,'position').*[1 1 0 1]+[0 0 14+tmp*49 0])
%
% color bars (just comment out next both lines)
% if tmp<0.5 ; clr = [tmp*2 1 0] ; else ; clr = [1 (1-tmp)*2 0] ; end% if
% set(h,'backgroundcolor',clr,'foregroundcolor',[.2 .2 .2])