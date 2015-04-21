/* etvc_mex.c - ensemble time varying convolution implemented as a matlab
    mex file.

    This file is part of MATLAB nlid toolbox.

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
 
    To compile type in Matlab command window: mex filename
        ex: mex etvc_mex.c

*/

#include <string.h>
#include "mex.h"


/* function prototypes: */
void mexFunction ( 
     int nlhs,                /* Number of left hand side arguments */
     mxArray *plhs[],          /* Length nlhs array for pntrs to lhs args */
     int nrhs,                /* Number of right hand side arguments */
     const mxArray *prhs[]           /* Length nlhs array for pntrs to rhs args */
   );

void time_varying_ensemble_convolution(
     double *xr1,                                  /* input ensemble */
     double *xr2,                              /* convolution kernel */
     double *xr3,                                 /* output ensemble */
     int num_realizations,        /* number of ensemble realizations */
     int length,                    /* number of points per ensemble */
     int n1,                                   /* lower filter bound */
     int n2,                                   /* upper filter bound */
     double dt                          /* the intersample increment */
   );

/* matlab user function input arguments */
#define X    prhs[0]
#define H    prhs[1]
#define DT   prhs[2]
#define N1   prhs[3]
#define N2   prhs[4]

/* matlab user function output arguments */
#define Y plhs[0]

/* matlab calling syntax:
   y = etvc_mex(x,h,dt,n1.n2);
   y - ensemble output, x ensemble convolved with timevarying kernel h
   x - the ensemble input
   h - the timevarying convolution kernel
   dt - the time increment (scalar)
   n1 - the lower convolution bound
   n2 - the upper convolution bound
*/

void mexFunction( int nlhs,           /* Number of left hand side arguments */
          mxArray *plhs[],       /* Length nlhs array for pntrs to lhs args */
          int nrhs,           /* Number of right hand side arguments */
          const mxArray *prhs[]        /* Length nlhs array for pntrs to rhs args */
        )
{

 


  /* first check the number of input/output arguments */
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments.");

  if (nrhs == 0)
      mexErrMsgTxt("Help message goes here");;
  if (nrhs != 5)
    mexErrMsgTxt("Incorrect number of input arguments.");

  /* check that the matrix dimensions agree */
  if ( mxGetM(X) != mxGetM(H))
    mexErrMsgTxt("The number of FIR's does not match the number of points in the ensemble inputs.");

  /* allocate space for the output matrix */
  Y = mxCreateDoubleMatrix (mxGetM(X), mxGetN(X), mxREAL);

  time_varying_ensemble_convolution( mxGetPr(X), mxGetPr(H), mxGetPr(Y),
			(int) mxGetN(X), 
            (int) mxGetM(X), (int)(*mxGetPr(N1)), 
            (int) *mxGetPr(N2), *mxGetPr(DT));

}  /* end of user_fnct() */
 

void time_varying_ensemble_convolution(
     double *x,                                    /* input ensemble */
     double *h,                                /* convolution kernel */
     double *y,                                   /* output ensemble */
     int num_realizations,        /* number of ensemble realizations */
     int length,                    /* number of points per ensemble */
     int tow1,                                 /* lower filter bound */
     int tow2,                                 /* upper filter bound */
     double dt                          /* the intersample increment */
     )

{
   int i,j,k,tow;
   double sum;
  
   for (k = 0; k <= num_realizations-1; k++) {
     for (i = 0; i <= length-1; i++) {
       sum = 0.0;
       /* An if is used to eliminate convolution "wrap around".  The
          input element is selected using tow (i-tow will always be
          between 0 and length-1), and the matching filter value is 
          indexed with j */ 
       for (tow = tow1, j = 0; tow <= tow2; tow++, j++) {
         if (i-tow >= 0  &&  i-tow <= length-1)
           sum += *(h + i + j*length) * 
                  *(x + (i-tow) + k*length); 
       }
       *(y + i + k*length) = sum*dt;
     }
   }

} /* end of function time_varying_ensemble_convolution() */

/* end of etvc_mex_v2.c */






