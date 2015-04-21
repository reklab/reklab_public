/*

  phix1y.c    mex-file to compute first-order cross-correlation in the 
              time domain
		 

  Usage:       phi = phix1y(x,y,hlen);


  David Westwick   November 16, 1993
 
*/


#include <math.h>
#include "mex.h"

/* Input Arguments */
#define	X_IN	prhs[0]
#define	Y_IN	prhs[1]
#define H_IN    prhs[2]


/* Output Arguments */
#define	PHI_OUT	plhs[0]


/*  macros */
#define	max(A, B)	((A) > (B) ? (A) : (B))
#define	min(A, B)	((A) < (B) ? (A) : (B))






/*  Function Prototypes  */

/*
    This routine computes the mean of a vector, x, and removes it.
    data_length is the length of the vector x

*/

void extern subtract_mean (
      double *x,				  /*  input channel  */
      unsigned int  data_length			  /* length of chanel */
    );







/*
    This is the main computation routine.   It computes the first-order
    cross-correlation between x and y, out to a memory length of lags.
    data_length is the length of the vectors.  The result is placed in
    a lags by lags matrix, phi.

    phi (i) = E [x(n-i) y(n)]
*/

void phixy(
	     double	         phi[],
	     double	         x[],
	     double              y[],
	     int                 lags, 
	     unsigned int        data_length
	)

              
{
  int i,j,k,lag_min;
  double sum,*ptr_i,*ptr_j,*ptr_k;
  
  
  for (i=0;i<lags;i++) {
	sum = 0.0;
	ptr_k = y+i;
	ptr_j = x;
	for (k=i;k<data_length;k++) {
	  sum += *ptr_k++ * *ptr_j++ ;
	}
	*(phi+i) = sum / data_length;
      }
}   /*   end of function phixy    */







/*
    This is the entry point for MATLAB

*/

#ifdef __STDC__
void mexFunction(
		 int		nlhs,
		 Matrix	        *plhs[],
		 int		nrhs,
		 Matrix	        *prhs[]
		 )
#else
mexFunction(nlhs, plhs, nrhs, prhs)
int nlhs, nrhs;
Matrix *plhs[], *prhs[];
#endif

{
  double           *phi;
  double           *x,*y,*hlen;
  double           *x0, *y0, *ptr_x, *ptr_y;
  unsigned int     lags, data_length, m, n;
 

  /* Check for proper number of arguments */
  
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phixy.m \n");
    printf("\n Usage:  phi = phixy(x,y,hlen); \n");
    printf("\t x:  input signal \n");
    printf("\t y:  output signal \n");
    printf("\t hlen:  memory length of phi \n");
  }
  else if (nrhs != 3) {
    mexErrMsgTxt("phix1y requires three input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("phix1y requires one output argument.");
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

    PHI_OUT = mxCreateFull(lags, 1, REAL);
    phi = mxGetPr(PHI_OUT);

  /*  Create and load work space arrays  */

    x0 = (double *) mxCalloc (2*data_length, sizeof(double));
    y0 = x0 + data_length;

    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));



  /* Do the actual computations in  other subroutines */

    subtract_mean (x0,data_length);
    subtract_mean (y0,data_length);
    phixy (phi,x0,y0,lags,data_length);
    return;
  }

}




















