/*
  cor_x3y.c    mex-file to compute third-order cross-correlation.

  Usage:       phi = corx3y(x,y,hlen);
  David Westwick   Original Version    August   24, 1993
                   Matlab 5 Port       December 14, 1999
                   Faster Computations December 17, 1999

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
*/

#include <math.h>
#include "mex.h"
/* Modified by Diego Guarin 2012/02/17 */  
#include <string.h>

/* Input Arguments */

#define	X_IN	prhs[0]
#define	Y_IN	prhs[1]
#define H_IN    prhs[2]

/* Output Arguments */

#define	PHI_OUT	plhs[0]

/*  macros Commented out Diego Guarin 2012/02/17 */

#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))



/*

    This is the main computation routine.   It computes the third-order
    cross-correlation between x and y, out to a memory length of lags.
    data_length is the length of the vectors.  The result is placed in
    a lags by lags by lags matrix, phi.

    phi (i,j,k) = E [x(n-i) x(n-j) x(n-k) y(n)]

*/


void phix3y(phi,x,y,lags,data_length)

	     double	         phi[];
	     double	         x[];
	     double          y[];
	     mwSize lags, data_length;
{

  mwSize i,j,k,n,packed_length,lags2;

  double temp1,temp2,*ptr_x,*ptr_x1,*ptr_x2,*ptr_x3;
  double *ptr_y,*packed,*phi1;
  packed_length = (lags * (lags+1) * (lags+2))/6;
  packed = (double *) mxCalloc (packed_length, sizeof(double));  
  ptr_y = y;
  ptr_x = x;

  /* note, ptr_x1 et all will occaisionally point to the lags elements 
     before the start of the x vector -- these have been initialized to 
*/
  for(n=0;n<data_length;n++,ptr_y++,ptr_x++){
    ptr_x1 = ptr_x;
    phi1 = packed;  
    for(i=0;i<lags;i++){
      ptr_x2 = ptr_x1;
      temp1 = *ptr_y * *ptr_x1--;
      for(j=i;j<lags;j++){
        ptr_x3 = ptr_x2;
        temp2 = temp1 * *ptr_x2--;
        for(k=j;k<lags;k++,phi1++)
          *phi1 += *ptr_x3-- * temp2;
      }
    }
  }

  phi1 = packed;
  lags2 = lags*lags;
  for (i=0;i<lags;i++) {
    for (j=i;j<lags;j++) {
      for (k=j;k<lags;k++,phi1++) {
        temp1 = *phi1/data_length;
        *(phi + lags2*i + lags*j + k) = temp1;
        *(phi + lags2*i + lags*k + j) = temp1;
        *(phi + lags2*j + lags*i + k) = temp1;
        *(phi + lags2*j + lags*k + i) = temp1;
        *(phi + lags2*k + lags*j + i) = temp1;
        *(phi + lags2*k + lags*i + j) = temp1;
      }
    }
  }
}                /*end of function phix3y */


/*    This is the entry point for MATLAB */
void mexFunction(nlhs, plhs, nrhs, prhs)
	int nlhs, nrhs;
	mxArray *plhs[];
	const mxArray *prhs[];

	{
  double           *phi;
  double           *x,*y,*hlen;
  double           *x0, *y0, *ptr_x, *ptr_y;
  mwSize     lags, data_length, m, n;

  /* Check for proper number of arguments */
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phix3y.m \n");
    printf("\n Usage:  phi = phix3y(x,y,hlen); \n");
    printf("\t x:  input signal \n");
    printf("\t y:  output signal \n");
    printf("\t hlen:  memory length of phi \n");

  }

  else if (nrhs != 3) {
    mexErrMsgTxt("phix3y requires three input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("phix3y requires one output argument.");
  } else {   /* proceed with operator */    

  /* Assign pointers to the various parameters */
    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    hlen = mxGetPr(H_IN);
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    data_length = max(m,n);
    lags = (unsigned int)*hlen;         /*Diego Guarin 2012/02*/
  /* Create a matrix for the return argument */

    PHI_OUT = mxCreateDoubleMatrix(lags*lags*lags, 1, mxREAL);
    phi = mxGetPr(PHI_OUT);
  /*  Create and load work space arrays  */

    /* the idea is to create a pad of lags zeros behind the start of x0, 
       so that we can avoid the array boundary checking in the main 
       computational routine */

    x0 = (double *) mxCalloc (2*data_length+lags, sizeof(double));
    x0 = x0+lags;
    y0 = x0 + data_length;
    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));

  /* Do the actual computations in  other subroutines */
    phix3y (phi,x0,y0,lags,data_length);
    return;
  }
}

