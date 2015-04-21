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
cd nlid_mex/src

mex -c  -largeArrayDims matrix_ref.c
mex  -largeArrayDims corx2y.c

% check for legacy OS that requires 8.3 filenames
if strcmpi(computer,'pcwin')
  mex corx3y.c matrix_ref.obj
  mex fast_pmpr.c matrix_ref.obj
elseif strcmpi(computer,'pcwin64')
  mex corx3y.c matrix_ref.obj
  mex fast_pmpr.c matrix_ref.obj
else
  mex  -largeArrayDims corx3y.c matrix_ref.o
  mex  -largeArrayDims fast_pmpr.c matrix_ref.o
end

switch lower(computer)
  case 'mac'
    !cp corx2y.mexmac ../../@cor/private
    !cp corx3y.mexmac ../../@cor/private
    !cp fast_pmpr.mexmac ../../@vseries/private
  case 'maci64' 
    !cp corx2y.mexmaci64 ../../@cor/private
    !cp corx3y.mexmaci64 ../../@cor/private
    !cp fast_pmpr.mexmaci64 ../../@vseries/private
  case 'glnx86'
    !cp corx2y.mexglx ../../@cor/private
    !cp corx3y.mexglx ../../@cor/private
    !cp fast_pmpr.mexglx ../../@vseries/private
    case 'lnx86'
    !cp corx2y.mexlx ../../@cor/private
    !cp corx3y.mexlx ../../@cor/private
    !cp fast_pmpr.mexlx ../../@vseries/private
  case 'sol2'
    !cp corx2y.mexsol ../../@cor/private
    !cp corx3y.mexsol ../../@cor/private
    !cp fast_pmpr.mexsol ../../@vseries/private
  case 'hpux'
    !cp corx2y.mexhpux ../../@cor/private
    !cp corx3y.mexhpux ../../@cor/private
    !cp fast_pmpr.mexhpux ../../@vseries/private
  case 'hp700'
    !cp corx2y.mexhp7 ../../@cor/private
    !cp corx3y.mexhp7 ../../@cor/private
    !cp fast_pmpr.mexhp7 ../../@vseries/private
  case 'ibm_rs'
    !cp corx2y.mexrs6 ../../@cor/private
    !cp orx3y.mexrs6 ../../@cor/private
    !cp fast_pmpr.mexrs6 ../../@vseries/private
  case 'sgi'
    !cp corx2y.mexsg ../../@cor/private
    !cp corx3y.mexsg ../../@cor/private
    !cp fast_pmpr.mexsg ../../@vseries/private
  case 'alpha'
    !cp corx2y.mexaxp ../../@cor/private
    !cp corx3y.mexaxp ../../@cor/private
    !cp fast_pmpr.mexaxp ../../@vseries/private 
  case 'pcwin'
    !copy corx2y.dll ..\..\@cor\private
    !copy corx3y.dll ..\..\@cor\private
    !copy fast_pmpr.dll ..\..\@vseries\private 
    case 'pcwin64'                                  %D. Guarin
    !copy corx2y.mexw64 ..\..\@cor\private       
    !copy corx3y.mexw64 ..\..\@cor\private
    !copy fast_pmpr.mexw64 ..\..\@vseries\private 
  otherwise    
    disp('unrecognized OS')
    disp('Move MEX files manually');
end

eval(['cd ' old]);
  



