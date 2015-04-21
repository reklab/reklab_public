function poly_gui(option,mdls,vafs)
% GUI application to choose the order of a polynomial
% 
% 
%  'option' is a control string.  It can take on the following values:
%    accept :  set the "done flag" and exit.
%    down   :  Decrease order by 1
%    init   :  Initialize the GUI, open figures, etc...
%    up     :  Increase order by 1
%

% Copyright 2003, Robert E Kearney and David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see ../../copying.txt and ../../gpl.txt 


% globals common @wbose/nlident
global  POLY_DONE POLY_ORDER

% globals common to previous calls to pdm_gui
global PLOY_VAFS POLY_MDLS POLY_FIGURE

switch option
  case 'accept'
    clf
    poly_gui('plot');
    POLY_DONE = 1;
    
  case 'init'
     POLY_DONE = 0;
     POLY_FIGURE = figure;
     POLY_MDLS = mdls;
     POLY_VAFS = vafs;
     [val,pos] = min(mdls);
     POLY_ORDER = pos - 1;
     poly_gui('plot');
     poly_gui('ui_controls');

   case 'plot'
     % plot the singular values, hiliting the current order.
     % also show the corresponding first and second-order kernels
     figure(POLY_FIGURE);
     orders = [1:length(POLY_MDLS)]'-1;
     
     subplot(211)
     plot(orders,POLY_VAFS,'*');
     title('Variance Accounted For');
     xlabel('Polynomial Order');
     hold on
     plot(POLY_ORDER,POLY_VAFS(POLY_ORDER+1),'m*');
     hold off     
     
     
     subplot(212)
     semilogy(orders,POLY_MDLS,'*');
     title('MDL Cost Function');
     xlabel('Polynomial Order');
     hold on
     plot(POLY_ORDER,POLY_MDLS(POLY_ORDER+1),'m*');
     hold off     
     
   case 'ui_controls'
     %
     %  Set up the UI controls
     %
     figure(POLY_FIGURE)
     subplot(211)
     accept_button = uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.75 0.65 0.12 0.0375], 'string',...
        'Accept','callback', 'poly_gui(''accept'');' );
     down_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.75 0.6 0.05 0.0375], 'string',...
        '<','callback', 'poly_gui(''down'');' );
      up_button =  uicontrol('style','pushbutton', 'units', 'normalized',...
        'pos', [0.82 0.6 0.05 0.0375], 'string',...
        '>','callback', 'poly_gui(''up'');' );
    
  case 'up'
    max_order = length(POLY_MDLS)-1;
    if POLY_ORDER < max_order
      POLY_ORDER = POLY_ORDER + 1;
      poly_gui('plot');
    end
  case 'down'
    if POLY_ORDER > 0
      POLY_ORDER = POLY_ORDER - 1;
      poly_gui('plot');
    end
    
  otherwise
    error('unrecognized option');
end

  


