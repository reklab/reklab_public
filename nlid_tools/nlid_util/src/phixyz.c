/*

  phixyz.c    mex-file to compute second-order cross-correlation.
		 

  Usage:       phi = phixyz(x,y,z,hlen);


  David Westwick  Version for DEC    August  24, 1993
                  Version for SUN    October 10, 1994
  REK Modifications for V5 and windows. 
 
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






/*
    This is the main computation routine.   It computes the second-order
    cross-cross-correlation between inputs x and y, and output z, out to a
    memory length of lags.
    "data_length" is the length of the vectors.  The result is placed in
    a lags by lags matrix, phi.

    phi (i,j) = E [x(n-i) y(n-j) z(n)]
*/

#ifdef __STDC__
void phixyz(
	     double	         phi[],
	     double	         x[],
	     double              y[],
	     double              z[],
	     int                 lags, 
	     unsigned int        data_length
	)
#else
phixyz(phi,x,y,z,lags,data_length)
	     double	         phi[];
	     double	         x[];
	     double              y[];
	     double              z[];
	     int                 lags; 
	     unsigned int        data_length;
#endif

              
{
  int i,j,k;
  double sum,*ptr_i,*ptr_j,*ptr_k;
  
  
  for (i=0;i<lags;i++) {
    for (j=0;j<lags;j++) {
      sum = 0.0;
      if ( i>j ) {
	ptr_k = z + i;
	ptr_j = y +i -j;
	ptr_i = x;
	for (k=i;k<data_length;k++) {
	  sum += *ptr_k++ * *ptr_i++ * *ptr_j++;
	}
      }
      else {
	ptr_k = z + j;
	ptr_j = y;
	ptr_i = x + j - i;
	for (k=j;k<data_length;k++) {
	  sum += *ptr_k++ * *ptr_i++ * *ptr_j++;
	}
      }	    
      *(phi+i+j*lags) = sum / data_length;
    }
  }
}   /*   end of function phixyz    */







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
  mxArray           *phi;
  double          *x,*y,*z,*hlen;
  double           *x0, *y0, *z0;
  double           *ptr_x, *ptr_y, *ptr_z;
  unsigned int     lags, data_length, m, n,i;
 

  /* Check for proper number of arguments */
  
  if (nrhs == 0) {     /*  print help message */
    printf("\n MEX-file implementation of phixyz.m \n");
    printf("\n Usage:  phi = phixyz(x,y,z,hlen); \n");
    printf("\t x:  first input signal \n");
    printf("\t y:  second input signal \n");
    printf("\t z:  output signal \n");
    printf("\t hlen:  memory length of phi \n");
  }
  else if (nrhs != 4) {
    mexErrMsgTxt("phixyz requires four input arguments.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("phixyz requires one output argument.");
  } else {   /* proceed with operator */    

  /* Assign pointers to the various parameters */
    printf("\nReady 1;\n");

    x = mxGetPr(X_IN);
    y = mxGetPr(Y_IN);
    z = mxGetPr(Z_IN);
    hlen = mxGetPr(H_IN);
  
    m = mxGetM(Y_IN);
    n = mxGetN(Y_IN);
    data_length = max(m,n);
    lags = (unsigned int)*hlen;

  /* Create a matrix for the return argument */
    printf("\nReady 2;\n");

    PHI_OUT = mxCreateDoubleMatrix(lags, lags, mxREAL);
    phi = mxGetPr(PHI_OUT);

  /*  Create and load work space arrays  */
    printf("\n 3 ;\n");

    x0 = (double *) mxCalloc (3*data_length, sizeof(double));
    y0 = x0 + data_length;
    z0 = y0 + data_length;

    printf("\n4;\n");

    ptr_x = (double *)memcpy(x0,x,data_length*sizeof(double));
    ptr_y = (double *)memcpy(y0,y,data_length*sizeof(double));
    ptr_z = (double *)memcpy(z0,z,data_length*sizeof(double));


    printf("\n6;\n");

  /* Do the actual computations in  other subroutines */


    printf("\n7;\n");
    subtract_mean (x0,data_length);
    subtract_mean (y0,data_length);
    subtract_mean (z0,data_length);

    printf("\nReady to call phixyz;\n");
    
    phixyz (phi,x0,y0,z0,lags,data_length);
    printf("\nphixyz Done ;\n");
    return;
  }
}
