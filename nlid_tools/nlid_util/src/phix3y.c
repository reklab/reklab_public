/*

  phi_x3y.c    mex-file to compute third-order cross-correlation.
		 

  Usage:       phi = phi_x3y(x,y,hlen);


  David Westwick   August 24, 1993
REK 20 May 97 Modifiied for matlab5 
REK 24 Sep 97 more changes
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


/*  function prototypes */

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

#ifdef __STDC__   
int extern matrix_ref (
		       int n,		              /* size of matrix */
		       int i,		  /* indeces of matrix  element */
		       int j,
		       int k
		       );
#else
int extern matrix_ref ();
		       int n;		              /* size of matrix */
		       int i;		  /* indeces of matrix  element */
		       int j;
		       int k;
#endif


/* end of function prototypes */



/*
    This is the main computation routine.   It computes the third-order
    cross-correlation between x and y, out to a memory length of lags.
    data_length is the length of the vectors.  The result is placed in
    a lags by lags by lags matrix, phi.

    phi (i,j,k) = E [x(n-i) x(n-j) x(n-k) y(n)]
*/

#ifdef __STDC__
void phix3y(
	     double	         phi[],
	     double	         x[],
	     double              y[],
	     int                 lags, 
	     unsigned int        data_length
	     )
#else
void phix3y(phi,x,y,lags,data_length)
	     double	         phi[];
	     double	         x[];
	     double              y[];
	     int                 lags; 
	     unsigned int        data_length;
#endif	     



{
  int i,j,k,n,lag_min;
  double sum;
  
  
  for (i=0;i<lags;i++) {
    for (j=0;j<= i;j++) {
      for (k=0;k <= j;k++) {
	sum = 0;
	for  (n=i;n<data_length;n++){
	  sum += *(y+n) * *(x+n-i) * *(x+n-j) * *(x+n-k);
	}
	sum = sum / data_length;
	*(phi + matrix_ref (lags,i,j,k)) = sum;
	*(phi + matrix_ref (lags,i,k,j)) = sum;
	*(phi + matrix_ref (lags,j,i,k)) = sum;
	*(phi + matrix_ref (lags,j,k,i)) = sum;
	*(phi + matrix_ref (lags,k,i,j)) = sum;
	*(phi + matrix_ref (lags,k,j,i)) = sum;
      }
    }
  }
}  /* end of function phix3y */







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
    lags = (unsigned int)*hlen;
 
  /* Create a matrix for the return argument */

    PHI_OUT = mxCreateDoubleMatrix(lags*lags*lags, 1, mxREAL);
    phi = mxGetPr(PHI_OUT);
  /*  Create and load work space arrays  */

    x0 = (double *) mxCalloc (2*data_length, sizeof(double));
    y0 = x0 + data_length;

    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));



  /* Do the actual computations in  other subroutines */

    subtract_mean (x0,data_length);
    subtract_mean (y0,data_length);
    phix3y (phi,x0,y0,lags,data_length);
    return;
  }
}
