function pinv_gui(option)
%  GUI application to choose deconvolution order in IRF estimation
%
%  'option' is a control string.  It can take on the following values:
%    accept :  set the "done flag" and exit.
%    down   :  Decrease order by 1
%    init   :  Initialize the GUI, open figures, etc...
%    up     :  Increase order by 1
%
%  pinv_gui expects to find global variables:
%    dtw_pinv_mdls   :  Minimum Descriptor Length cost function for all
%                       orders.
%    dtw_pinv_IRFs   :  Matrix of IRFs where the n'th column contains the
%                       IRF ontained using n singular vectors.  Initially,
%                       this matrix will contain all NaNs, except for the
%                       column corresponding to the minimum MDL.  As IRFs
%                       are needed for display, they will be computed.
%    dtw_pinv_Bounds :  Matrix of error bounds. Each column contains the
%                       confidence bound associated with the corresponding
%                       IRF in dtw_pinv_IRFs.  Initially, contains all NaNs,
%                       except for the minimum MDL....
%
%  pinv_gui returns the following globals
%    dtw_pinv_done   :  flag to indicate that the choice has been made
%    dtw_pinv_order  :  number of singular vectors included in final
%                       IRF.
%
% Copyright 1994-2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 
global dtw_mdls dtw_pinv_IRFs dtw_pinv_Bounds
global dtw_pinv_done dtw_pinv_order dtw_pinv_numlags
global dtw_pinv_irf dtw_pinv_bound
global dtw_pinv_figure dtw_pinv_marker
if strcmp(option,'accept')
  clf
  pinv_gui('plot');
  dtw_pinv_done = 1;
  return
  
elseif strcmp(option,'down')
  if dtw_pinv_order > 1
    dtw_pinv_order = dtw_pinv_order - 1;
    pinv_gui('plot');
  end
  return
  
elseif strcmp(option,'down5')
  if dtw_pinv_order > 5
    dtw_pinv_order = dtw_pinv_order - 5;
    pinv_gui('plot');
  end
  return
  
  
elseif strcmp(option,'init')
  dtw_pinv_done = 0;
  dtw_pinv_IRFs = nan*zeros(dtw_pinv_numlags);
  dtw_pinv_Bounds =  dtw_pinv_IRFs;
  if isempty(dtw_pinv_figure)
    dtw_pinv_figure = figure;
  end
  pinv_gui('plot');
  pinv_gui('ui_controls');
  return
  
elseif strcmp(option,'plot')
  %
  %  Plot graphs of the MDL and the automatically chosen IRF
  %
  mdls = dtw_mdls;
  %
  %  Check to see if the appropriate column has been generated
  %
  if isnan(dtw_pinv_IRFs(1,dtw_pinv_order))
    %disp(['computing IRF of order ' num2str(dtw_pinv_order)])
    compute_pinv('irf');
    dtw_pinv_IRFs(:,dtw_pinv_order) = dtw_pinv_irf;
    dtw_pinv_Bounds(:,dtw_pinv_order) = dtw_pinv_bound;
  else
    dtw_pinv_irf = dtw_pinv_IRFs(:,dtw_pinv_order);
    dtw_pinv_bound = dtw_pinv_Bounds(:,dtw_pinv_order);
  end
  irf = dtw_pinv_irf;
  figure(dtw_pinv_figure);
  subplot(211)
  semilogy(mdls);
  title('Minimim Descriptor Length Cost Function');
  xlabel('Pseudo-Inverse Order');
  hold on
  dtw_pinv_marker = plot(dtw_pinv_order,mdls(dtw_pinv_order),'m*');
  hold off
  subplot(212);
  plot(irf);
  if isempty(dtw_pinv_Bounds)
    title('IRF estimate');
  else
    bound = dtw_pinv_bound;
    hold on
    plot(irf + bound,'m');
    plot(irf - bound,'m');
    hold off
    title('IRF estimate between confidence bounds');
  end
  return
  
elseif strcmp(option,'ui_controls')
  %
  %  Set up the UI controls
  %
  figure(dtw_pinv_figure)
  subplot(211)
  accept_button = uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.75 0.85 0.12 0.0375], 'string',...
        'Accept','callback', 'pinv_gui(''accept'');' );
  down_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.75 0.8 0.05 0.0375], 'string',...
        '<','callback', 'pinv_gui(''down'');' );
  up_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.82 0.8 0.05 0.0375], 'string',...
        '>','callback', 'pinv_gui(''up'');' );
  down5_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.75 0.75 0.05 0.0375], 'string',...
        '-5','callback', 'pinv_gui(''down5'');' );
  up5_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos',[0.82 0.75 0.05 0.0375], 'string',...
        '+5','callback', 'pinv_gui(''up5'');' );
    
  return
elseif strcmp(option,'up')  
  if dtw_pinv_order < length(dtw_pinv_IRFs(:,1));
    dtw_pinv_order = dtw_pinv_order + 1;
    pinv_gui('plot');
  end
  return
elseif strcmp(option,'up5')  
  if dtw_pinv_order < (length(dtw_pinv_IRFs(:,1))-4);
    dtw_pinv_order = dtw_pinv_order + 5;
    pinv_gui('plot');
  end
  return
end
