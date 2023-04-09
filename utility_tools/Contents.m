% UTILITY_TOOLS
%
% Files
%   any2csv               - any2text writes contents of s to a nicely formated csv-file
%   arg_help              - display help for arguments
%   arg_parse             - arg_parse ( options, varargin) -Parse and decode variable number of inputs
%   arg_parse_c           - arg_parse ( caseMode, options, varargin) -Parse and decode variable number of inputs with case control
%   arg_parseAgent        - arg_parse ( options, varargin) -Parse and decode variable number of inputs
%   chop                  - Round matrix elements.
%   CMB_time              - Time and date utilities
%   cseq2eseq.old         - catSignal2Events - convert a categorical sequence  to event sequence
%   cseqPlot              - PLOTPATTERNS Plots the sequence of events P in a cseq File 
%   deblankl              - Remove leading blanks.
%   del                   - This function adds a delay to an input signal (x) that is not a nldat object
%   delimn_read           - 
%   dockFig               - dockFig - dock/undock figure(1)
%   efm                   - efm - class for  manipulation of efm signals
%   eseqPlot              - PLOTPATTERNS Plots the sequence of events P
%   fgzw                  - compute 2-nd order IRF for k,z,w
%   fgzwp                 - compute 3rd order IRF for k,z,w,p
%   fig_mod_old           - function fig_mod  (h, 'option' ,'value' );
%   figMod                - function fig_mod  (h, 'option' ,'value' );
%   figPrint              - figPrint(figNum, fileName, printFlag, newFileFlag, pageOrient, formatType)
%   filter_ts             - performs one and two-sided filtering.
%   fkbi                  - compute second-order IRF for K,B,I
%   fkzw                  - compute 2-nd order IRF for k,z,w
%   flb2mat               - FLBTOMAT_M.M ...
%   flbio                 - - perform io operations on flb files
%   flbtomat              - read data from a  FILELB file
%   footer                - printer a footer message at the bottom of the screen
%   gen_pralds            - This function generates a desired number of realizations from a position/velocity
%   gen_pralds_Test       - 
%   guiInput_d            - guiInput_d  gui based interqactive input of numeric values with defaults and  range checking
%   guiInput_l            - input logical y/n answer
%   guiInput_List         - - select from a strong from a list
%   guiInput_s            - input a string from guiwith optional  default and list
%   hindex                - hindex (oper, ht, varargin) haded index array operations
%   html_token            - [token]=html_token ( S, start_string, stop_string, string_flag, ntoken)
%   idxfun                - idxfun -  support indexed arrays for clusters
%   input_d               - input_d  interqactive input of numeric values with defaults and  range checking
%   input_d1              - decode a  string with default values and option range checking
%   input_l               - input logical y/n answer
%   input_s               - input a string with optional  default and list
%   input_sc              - input a string with optional  default and list
%   inputsdlg             - Enhanced input dialog box supporting multiple data types
%   inputsdlg_demo        - INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)
%   install_nlid_mexfiles - compiles and installs mexfiles for NLID toolbox
%   invert_irf            - Inverts an impulse response function in the frequency domain
%   MakeContents          - makes Contents file in current working directory.
%   mat2flb               - ...
%   menu_reklab           - MENU   Generate a menu of choices for user input.
%   mm                    - Memory Monitor displays runtime memory information
%   multi_tcheb           - usage W = multi_herm2(V, max_order);
%   newfcn                - create a MATLAB function with entered filename
%   newPert               - This function is intended to generate the novel input perturbation signal (velocity perturbation) with a desired PDF for identification of
%   ossep                 - 
%   pmatch                - Case insensitive string matching with automatic completion.
%   pnvmain               - [PROPS,ASGNVALS] = pnv (SYS,'names', 'true')  returns the list PROPS of
%   printFigures          - printFigures (figureList, fileName)
%   pvdisp                - PVPDISP Displays properties PROPS and their values VALUES 
%   quote                 - 
%   rk_split              - split a string into components based on delimiter
%   rsync                 - - matlab wrapper for rsyn  
%   strcvrt               - 
%   streamer              - Titles for an entire figure.
%   struct2txt            - struct_table - output a structure to a csv file 
%   struct_disp           - struct_disp - display  struct with heading $Revision $
%   struct_print          - print table with heading
%   struct_table          - struct_table - output a structure to a html table
%   struct_table_old      - struct_table - output a structure to a html table
%   structop              - SS = StructOP (OPER,  S, varargin )
%   subplotLayout         - Select subplot layout for mulutiple plots
%   subset_fun            - subset_fun - supports subset operatins
%   table_print           - print table with heading
%   tv_2                  - Generates a time varying impusle response of second-order low-pass filters (IRF's). 
%   tv_2nd_order_LP       - Generates matrix of second-order low-pass filters (IRF's). 
%   tv_hammer_system      - hammer = hammer_system(duration,dt,verbose)
%   uigetdir2             - Standard dialog box for selecting a directory which remembers
%   uigetfile2            - Standard open file dialog box which remembers last opened folder
%   uiputfile2            - Standard save file dialog box which remembers last opened folder
%   unixt2datestr         - Convert a unix time stamp to a date time str
%   varg                  - examine varargin and correct for nesting calls
%   xAxisPanZoom          - xAxisPanZoom - enable linked pan and zoom for linked axes 
%   xticklabel_rotate     - hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
