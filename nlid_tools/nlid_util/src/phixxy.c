/*

  phixxy.c    mex-file to compute second-order cross-correlation.
		 

  Usage:       phi = phixxy(x,y,hlen);


  David Westwick  Version for DEC    August  23, 1993
                  Revised for SUN    October 10, 1994
                  REvised for Windoes 24 Sept 97 REK 
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

#ifdef __STDC__
void extern subtract_mean (
      double *x,				  /*  input channel  */
      unsigned int  data_length			  /* length of chanel */
    );
#else
extern void subtract_mean ();
      double *x;
      unsigned int data_length;
#endif






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
	     unsigned int        data_length
	)
#else
phixxy(phi,x,y,lags,data_length)
	     double	         phi[];
	     double	         x[];
	     double              y[];
	     int                 lags; 
	     unsigned int        data_length;
#endif


              
{
  int i,j,k,lag_min;
  double sum,*ptr_i,*ptr_j,*ptr_k;
  
  
  for (i=0;i<lags;i++) {
    for (j=0;j<lags;j++) {
      if ( i>j ) {
	*(phi+i+j*lags) = *(phi+j+i*lags);
      }   else   {	    
	sum = 0.0;
	ptr_k = y+j;
	ptr_j = x;
	ptr_i = x + j - i;
	for (k=j;k<data_length;k++) {
	  sum += *ptr_k++ * *ptr_i++ * *ptr_j++;
	}
	*(phi+i+j*lags) = sum / data_length;
      }
    }
  }
}   /*   end of function phixxy    */







/*
    This is the entry point for MATLAB

*/


void mexFunction(
		 int		nlhs,
		 mxArray	        *plhs[],
		 int		nrhs,
		 const mxArray	        *prhs[]
		 )

{
  double           *phi;
  double           *x,*y,*hlen;
  double           *x0, *y0, *ptr_x, *ptr_y;
  unsigned int     lags, data_length, m, n;
 

  /* Check for proper number of arguments */
  
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phixxy.m \n");
    printf("\n Usage:  phi = phixxy(x,y,hlen); \n");
    printf("\t x:  input signal \n");
    printf("\t y:  output signal \n");
    printf("\t hlen:  memory length of phi \n");
  }
  else if (nrhs != 3) {
    mexErrMsgTxt("phixxy requires three input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("phixxy requires one output argument.");
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

  /*  Create and load work space arrays  */

    x0 = (double *) mxCalloc (2*data_length, sizeof(double));
    y0 = x0 + data_length;

    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));



  /* Do the actual computations in  other subroutines */

    subtract_mean (x0,data_length);
    subtract_mean (y0,data_length);

    phixxy (phi,x0,y0,lags,data_length);
    return;
  }
}
