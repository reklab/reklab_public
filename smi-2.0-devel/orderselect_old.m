function order=orderselect(S,method)
% orderselect Gives an estimate of the order of a system based on the 
%             singular values. Two methods are available:
%             manual and largest-gap. The manual method shows a
%             semi-logarithmic plot of the singular values and lets the user
%             manually choose the order by inspection. The largest-gap
%             method simply selects the order to be the one after which the
%             largest gap occures in the semi-logarithmic plot.
%
% Syntax:
%             order=orderselect(S,method)
%
% Inputs:
%   S         Vector with singular values, obtained from one of the
%             preprocessor functions  dordxx(dordom, dordpi, etc).
%   method    Method for selecting the order. Possible methods are:
%             ``manual'' and ``largest-gap''. Default is ``manual''.
%
% Output: 
%   order     Selected model order.
%
% see also: moesp, dordpo, dordpi

if nargin==1
  order=orderselect(S,'manual');
end

if nargin==2
  if length(method)<3
    error('unknown method')
  end
  if length(method)>3
    method=method(1:3);
  end
  if method=='man'
    mygui=orderselect(S,'gui');
    busy_flag=1;
    while busy_flag==1
      pause(0.5);
      Data=get(mygui,'Userdata');
      busy_flag=Data.busy_flag;
    end
    order=Data.order;
    close(mygui);
    %fig=order;
  elseif method=='lar'
    % largest gap method
    dS=diff(log(S));
    order=find(dS==min(dS));
  elseif method=='gui'
    % initialize gui for manual method
    fig = figure(...
	'Tag','orderselectgui',...
	'Color',[0.8 0.8 0.8], ...
	'Units','Normalized');
    axis1=axes(...
	'Parent',fig, ...
	'Tag','orderselectfig',...
	'Units','Normalized', ...
	'ButtonDownFcn','orderselect([],''bdf'');', ...
	'Position',[0.05 0.15 0.93 0.8]);
    buttonleft=uicontrol(...
	'Parent',fig, ...
	'Tag','buttonleft',...
	'Style','pushbutton', ...
	'Units','Normalized', ...
	'Position',[0.05 0.02, 0.16 0.08], ...
	'String','<<', ...
	'Callback','orderselect(1,''callback'');');
    buttonrightt=uicontrol(...
	'Parent',fig, ...
	'Style','pushbutton', ...
	'Units','Normalized', ...
	'Position',[0.25 0.02, 0.16 0.08], ...
	'String','>>', ...
	'Callback','orderselect(2,''callback'');');
    buttonok=uicontrol(...
	'Parent',fig, ...
	'Style','pushbutton', ...
	'Units','Normalized', ...
	'Position',[0.62 0.02, 0.16 0.08], ...
	'String','OK', ...
	'Callback','orderselect(3,''callback'');');
    buttoncnl=uicontrol(...
	'Parent',fig, ...
	'Style','pushbutton', ...
	'Units','Normalized', ...
	'Position',[0.82 0.02, 0.16 0.08], ...
	'String','Cancel', ...
	'Callback','orderselect(4,''callback'');');
    title1 = text('Parent',axis1, ...
	'Units','Normalized',...
	'Color',[0 0 0], ...
	'HandleVisibility','off', ...
	'HorizontalAlignment','center', ...
	'Position',[0.05 0.04 0.7],...
	'VerticalAlignment','cap');
    
    % Creat userdata
    dS=diff(log(S(1:end-1))); %forget last value of S, it is often bogus
    order=find(dS==min(dS));
    Data=struct('busy_flag',1,'S',S,'order',order);
    set(fig,'Userdata',Data); 
    % plot data
    hold off
    H=semilogy(S,'x');
    set(H,'ButtonDownFcn','orderselect(2,''bdf'');');
    hold on
    plot(order,S(order),'r.');
    Fig=get(fig,'Children');Fig=Fig(5);
    Dot=get(Fig,'Children');Dot=Dot(1);
    oldMarkersize= get(Dot,'Markersize');
    set(Dot,'Markersize',oldMarkersize*4);  
    title(['Selected order: ',num2str(order)]);  
    %return figure handle
    order=fig; 
  elseif method=='cal' %callback
    cbno=S;
    fig=gcbf;
    Data=get(fig,'Userdata');
    S=Data.S;
    order=Data.order;
    maxorder=length(S);
    if (cbno==1)
      %callback of buttonleft
       order=order-1;
      if order==0;
	order=1 ;
      end
      hold off
      H=semilogy(S,'x');
      set(H,'ButtonDownFcn','orderselect(2,''bdf'');');
      hold on
      plot(order,S(order),'r.');
      Fig=get(fig,'Children');Fig=Fig(5);
      Dot=get(Fig,'Children');Dot=Dot(1); 
      oldMarkersize= get(Dot,'Markersize');
      title(['Selected order: ',num2str(order)]);
      set(Dot,'Markersize',oldMarkersize*4);  
      Data.order=order; 
      set(fig,'Userdata',Data); 
    elseif (cbno==2)
      %callback of buttonright
      order=order+1;
      if order>maxorder;
	order=maxorder ;
      end
      hold off
      H=semilogy(S,'x');
      set(H,'ButtonDownFcn','orderselect(2,''bdf'');');
      hold on
      plot(order,S(order),'r.');
      Fig=get(fig,'Children');Fig=Fig(5);
      Dot=get(Fig,'Children');Dot=Dot(1); 
      oldMarkersize= get(Dot,'Markersize');
      set(Dot,'Markersize',oldMarkersize*4);  
      title(['Selected order: ',num2str(order)]);
      Data.order=order;
      set(fig,'Userdata',Data); 
    elseif (cbno==3)
      %callback of buttonok
      Data.busy_flag=0;
      set(fig,'Userdata',Data); 
    elseif (cbno==4)
      % callback of buttoncnl
      Data.order=[];
      Data.busy_flag=0;
      set(fig,'Userdata',Data); 
    end %cbno
  elseif method=='bdf'
    % buttondown function of figure
    fig=gcbf;
    Data=get(fig,'Userdata');
    S=Data.S;
    order=Data.order;
    Currpoint=get(gca,'CurrentPoint');
    order=round(Currpoint(1,1));
    hold off
    H=semilogy(S,'x');
    set(H,'ButtonDownFcn','orderselect(2,''bdf'');');
    hold on;
    plot(order,S(order),'r.');
    Fig=get(fig,'Children');Fig=Fig(5);
    Dot=get(Fig,'Children');Dot=Dot(1); 
    oldMarkersize= get(Dot,'Markersize');
    title(['Selected order: ',num2str(order)]);
    set(Dot,'Markersize',oldMarkersize*4);  
    title(['Selected order: ',num2str(order)]);
    Data=get(fig,'Userdata');
    Data.order=order;
    set(fig,'Userdata',Data); 
  end %method
end %nargin






