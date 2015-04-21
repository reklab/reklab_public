%
% Build all NLID mex files
%


%return
%disp('Make k3_sym')
%mex k3_sym.c matrix_ref.c
%disp('make k3_gen');
%mex k3_gen.c matrix_ref.c
%disp('Make k32_gen')
%mex k32_gen.c matrix_ref.c
%disp('Make gconv31')
%mex gconv31.c matrix_ref.c
%disp ('Make mult32')
%mex mult32.c matrix_ref.c
%disp ('Make mult3_101')
%mex mult3_101.c matrix_ref.c 
disp('making mex files') 
mex -c  -largeArrayDims matrix_ref.c
mex  -largeArrayDims corx2y.c
mex  -largeArrayDims corx3y.c matrix_ref.o
mex  -largeArrayDims fast_pmpr.c matrix_ref.o
mex phix2yz.c matrix_ref.c subtract_mean.c
mex phix3y.c matrix_ref.c subtract_mean.c
mex phixxy.c subtract_mean.c 
mex phixyz.c subtract_mean.c

option=input_l('Install NLID mex files ','y')
if option,
  echo on
  !delete ..\*.dll
  !move *.dll ..
  echo off
end

