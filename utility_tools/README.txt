This archive contains the MATLAB nonlinear system identification
(NLID) toolbox.  The algorithms implemented in this toolbox are 
described in
  Westwick, DT and Kearney, RE, "Identification of Nonlinear
  Physiological Systems: Theory and Practice", IEEE Engineering in
  Medicine and Biology Book Series, IEEE Press/John Wiley & Sons,
  2003.  ISBN: 0471274569.

LICENSE

The nlid toolbox is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

The nlid toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the nlid toolbox; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307  USA


CONTENTS

The archive should contain 3 subdirectories

nlid_book
  nlid_book/examples -- contains the code used to generate figures in
  Westwick and Kearney, 2003.

  nlid_book/problems -- contains solutions to all of the computer
  exercises in Westwick and Kearney, 2003.

nlid_tools -- contains the nlid toolbox routines.

utility_tools -- contains some support routines for the nlid toolbox.

INSTALLATION

Place the nlid_tools and utility_tools directories wherever MATLAB
toolboxes are stored on your system, and add both of them to your 
MATLAB search path. For example, if nlid_tools and utility_tools had
been placed in the directory /home/westwick/matlab/m_files, then
the commands:

path('/home/westwick/matlab/m_files/utility_tools',path);
path('/home/westwick/matlab/m_files/nlid_tools',path);

would add both these directories to the front of the search path (so
that they are searched before the rest of the path).  These commands 
may be added to your startup.m file, so that the NLID toolbox is
automatically enabled when MATLAB is started.

Several of the more computationally intense routines are implemented
as MEX files.  MEX files for Windows, Linux and MAC OSX have been
included in the archive.  Users of other operating systems will have
to build these files themselves.  

To build the MEX files, set up MATLAB to compile C/C++ MEX files.
Start MATLAB, and add the nlid_tools and utility_tools to the search
path.  Run the command install_nlid_mexfiles.m  The mex files have been
built using the GCC.  The behaviour of other complilers has not been
tested. 