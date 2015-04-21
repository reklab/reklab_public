/*

  phixyz.c    mex-file to compute third-order cross-cross-correlation.
		 

  Usage:       phi = phix2yz(x,y,z,hlen);


  David Westwick   August 24, 1993
  Revised for K+R compatibility (required for Sun SparcStations
  October 10, 1994 
  20 May 97 REK revised for matlab5
  24 Sep 97 REK for MSVC++
*/


#include <math.h>
#include "mex.h"

/* Input Arguments */
#define	X_IN	prhs[0]
#define	Y_IN	prhs[1]
#define Z_IN    prhs[2]
#define H_IN    prhs[3]


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
    cross-cross-correlation between inputs x and y, and output z, out to a
    memory length of lags.  The correlation is second-order in x and first
    order in y.
    "data_length" is the length of the vectors.  The result is placed in
    a lags by lags matrix, phi.

    phi (i,j,k) = E [x(n-i) x(n-j) y(n-k) z(n)]
*/

#ifdef __STDC__
void phix2yz(
	     double	         phi[],
	     double	         x[],
	     double              y[],
	     double              z[],
	     int                 lags, 
	     unsigned int        data_length
	)
#else
void phix2yz(phi,x,y,z,lags,data_length)
	     double	         phi[];
	     double	         x[];
	     double              y[];
	     double              z[];
	     int                 lags; 
	     unsigned int        data_length;
#endif


{
   int i,j,k,n,lag_min;
   double sum;
  
  
   for (i=0;i<lags;i++) {
     for (j=0;j<= i;j++) {
       for (k=0;k < lags;k++) {
	 lag_min = i;
	 if (k>i) {
	   lag_min = k;
	 }
	 sum = 0;
	 for  (n=lag_min;n<data_length;n++){
	   sum += *(z+n) * *(x+n-i) * *(x+n-j) * *(y+n-k);
	 }
	 sum = sum / data_length;
	 *(phi + matrix_ref (lags,i,j,k)) = sum;
	 *(phi + matrix_ref (lags,j,i,k)) = sum;
       }
     }
   }
 }  /* end of function corr_21 */







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
  double           *x,*y,*z,*hlen;
  double           *x0, *y0, *z0;
  double           *ptr_x, *ptr_y, *ptr_z;
  unsigned int     lags, data_length, m, n,i;
 

  /* Check for proper number of arguments */
  
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phix2yz.m \n");
    printf("\n Usage:  phi = phix2yz(x,y,z,hlen); \n");
    printf("\t x:  first input signal \n");
    printf("\t y:  second input signal \n");
    printf("\t z:  output signal \n");
    printf("\t hlen:  memory length of phi \n");
  }
  else if (nrhs != 4) {
    mexErrMsgTxt("phix2yz requires four input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("phix2yz requires one output argument.");
  } else {   /* proceed with operator */    

  /* Assign pointers to the various parameters */

    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    z = mxGetPr(Z_IN);
    hlen = mxGetPr(H_IN);
  
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    data_length = max(m,n);
    lags = (unsigned int)*hlen;

  /* Create a matrix for the return argument */

    PHI_OUT = mxCreateDoubleMatrix(lags*lags*lags, 1, mxREAL);
    phi = mxGetPr(PHI_OUT);

  /*  Create and load work space arrays  */

    x0 = (double *) mxCalloc (3*data_length, sizeof(double));
    y0 = x0 + data_length;
    z0 = y0 + data_length;

    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));
    ptr_z = (double *)memcpy(z0,z,data_length*sizeof(double));

  /* Do the actual computations in  other subroutines */

    subtract_mean (x0,data_length);
    subtract_mean (y0,data_length);
    subtract_mean (z0,data_length);
    phix2yz (phi,x0,y0,z0,lags,data_length);
    return;
  }
}


