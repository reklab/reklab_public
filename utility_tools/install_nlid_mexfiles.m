function install_nlid_mexfiles
% compiles and installs mexfiles for NLID toolbox


% Copyright 2003, Robert E Kearney, David T Westwick
% This file is part of the nlid toolbox, and is released under the GNU 
% General Public License For details, see copying.txt and gpl.txt 


% find where the nlid_tools are installed.
s = which('nlm');
index = find(s=='@');
base = s(1:index-1);

% store the current working directory
old = pwd;

% change directories to the nlid_tools root
eval(['cd ' base]);
cd nlid_util/src

mex -c matrix_ref.c
mex corx2y.c

% check for legacy OS that requires 8.3 filenames
if strcmp(lower(computer),'pcwin')
  mex corx3y.c matrix_ref.obj
  mex fast_pmpr.c matrix_ref.obj
else
  mex corx3y.c matrix_ref.o
  mex fast_pmpr.c matrix_ref.o
end

switch lower(computer)
  case 'mac'
    !mv corx2y.mexmac ../../@cor/private
    !mv corx3y.mexmac ../../@cor/private
    !mv fast_pmpr.mexmac ../../@vseries/private
  case 'glnx86'
    !mv corx2y.mexglx ../../@cor/private
    !mv corx3y.mexglx ../../@cor/private
    !mv fast_pmpr.mexglx ../../@vseries/private
    case 'lnx86'
    !mv corx2y.mexlx ../../@cor/private
    !mv corx3y.mexlx ../../@cor/private
    !mv fast_pmpr.mexlx ../../@vseries/private
  case 'sol2'
    !mv corx2y.mexsol ../../@cor/private
    !mv corx3y.mexsol ../../@cor/private
    !mv fast_pmpr.mexsol ../../@vseries/private
  case 'hpux'
    !mv corx2y.mexhpux ../../@cor/private
    !mv corx3y.mexhpux ../../@cor/private
    !mv fast_pmpr.mexhpux ../../@vseries/private
  case 'hp700'
    !mv corx2y.mexhp7 ../../@cor/private
    !mv corx3y.mexhp7 ../../@cor/private
    !mv fast_pmpr.mexhp7 ../../@vseries/private
  case 'ibm_rs'
    !mv corx2y.mexrs6 ../../@cor/private
    !mv corx3y.mexrs6 ../../@cor/private
    !mv fast_pmpr.mexrs6 ../../@vseries/private
  case 'sgi'
    !mv corx2y.mexsg ../../@cor/private
    !mv corx3y.mexsg ../../@cor/private
    !mv fast_pmpr.mexsg ../../@vseries/private
  case 'alpha'
    !mv corx2y.mexaxp ../../@cor/private
    !mv corx3y.mexaxp ../../@cor/private
    !mv fast_pmpr.mexaxp ../../@vseries/private 
  case 'pcwin'
    !move corx2y.dll ..\..\@cor\private
    !move corx3y.dll ..\..\@cor\private
    !move fast_pmpr.dll ..\..\@vseries\private 
  otherwise    
    disp('unrecognized OS')'
    disp('Move MEX files manually');
end

eval(['cd ' old]);
  



