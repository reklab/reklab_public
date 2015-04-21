/*

  corx2y.c    mex-file to compute second-order cross-correlation.
              without subtracting means (as for FOA)

  Usage:       phi = corx2y(x,y,hlen);


  David Westwick  Version for DEC    August  23, 1993
                  Revised for SUN    October 10, 1994
                  Port to Matlab 5   December 10, 1999
                  Fast Computations  December 13, 1999

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
#include "matrix.h"

/* Input Arguments */
#define	X_IN	prhs[0]
#define	Y_IN	prhs[1]
#define H_IN    prhs[2]


/* Output Arguments */
#define	PHI_OUT	plhs[0]


/*  macros  Commented out Diego Guarin 2012/02/17 */
#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))






/*  Function Prototypes  */

/*
    This routine computes the mean of a vector, x, and removes it.
    data_length is the length of the vector x

*/






/*
    This is the main computation routine.   It computes the second-order
    cross-correlation between x and y, out to a memory length of lags.
    data_length is the length of the vectors.  The result is placed in
    a lags by lags matrix, phi.

    phi (i,j) = E [x(n-i) x(n-j) y(n)]
*/

#ifdef __STDC__
void phixxy(
	     double	         phi[],
	     double	         x[],
	     double              y[],
	     int                 lags, 
	     int        data_length
	)
#else
phixxy(phi,x,y,lags,data_length)
	     double	         phi[];
	     double	         x[];
	     double              y[];
	     mwSize                 lags, data_length;
#endif

{
  mwSize i,j,k;
  double temp,*ptr_x,*ptr_y,*phi1;
  

  ptr_y = y;
  for(k=0;k<lags;k++,ptr_y++){
    for(i=0;i<=k;i++){
      ptr_x = x +k - i;
      temp = *ptr_y * *ptr_x;
      phi1 = phi + i*(lags+1);
      for(j=i;j<=k;j++,phi1++)
        *phi1 += temp * *ptr_x--;
    }
  }
  for(k=lags;k<data_length;k++,ptr_y++){
    phi1 = phi;
    for(i=0;i<lags;i++){
      ptr_x = x +k - i;
      temp = *ptr_y * *ptr_x;
      phi1 = phi1 + i;
      for(j=i;j<lags;j++,phi1++)
        *phi1 += temp * *ptr_x--;
    }
  }

  phi1 = phi;
  for (i=0;i<lags;i++) {
    *phi1 = *phi1/data_length;
    ptr_x = phi + i;
    ptr_y = phi + i*lags;
    for (j=0;j<i;j++){
      *ptr_x = *ptr_x/ data_length;
      *ptr_y++ = *ptr_x;
      ptr_x += lags;
    }
    phi1 += lags+1;
  }

}   /*   end of function phixxy    */







/*
    This is the entry point for MATLAB

*/

#ifdef __STDC__
void mexFunction(
		 int		nlhs,
		 mxArray	*plhs[],
		 int		nrhs,
		 const mxArray  *prhs[]
		 )
#else
void mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
mxArray *plhs[];
const mxArray *prhs[];
#endif

{
  double           *phi;
  double           *x,*y,*hlen;
  double           *x0, *y0, *ptr_x, *ptr_y;    
  mwSize   m,n, lags, data_length;
 

  /* Check for proper number of arguments */
  
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phixxy.m \n");
    printf("\n Usage:  phi = corxxy(x,y,hlen); \n");
    printf("\t x:  input signal \n");
    printf("\t y:  output signal \n");
    printf("\t hlen:  memory length of phi \n");
  }
  else if (nrhs != 3) {
    mexErrMsgTxt("corxxy requires three input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("corxxy requires one output argument.");
  } else {   /* proceed with operator */    

  /* Assign pointers to the various parameters */

    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    hlen = mxGetPr(H_IN);
  
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    data_length = max(m,n);
    lags = (unsigned int)*hlen;

  /* Create a matrix for the return argument */

    PHI_OUT = mxCreateDoubleMatrix(lags, lags, mxREAL);
    phi = mxGetPr(PHI_OUT);


  /* Do the actual computations in  other subroutines */

    phixxy (phi,x,y,lags,data_length);
    return;
  }
}
