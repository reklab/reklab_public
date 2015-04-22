/*



  fast_pmpr.c    mex-file to compute the pmpr matrix (X'X really)

              used by the fast orthogonal algorithm



  Usage:       pmpr = fast_pmpr(phiuy,phiu2y,phiu3y,u);



  David Westwick 25 Oct 1997





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



/* Input Arguments */

#define	P1_IN	prhs[0]

#define	P2_IN	prhs[1]

#define P3_IN   prhs[2]

#define	U_IN	prhs[3]





/* Output Arguments */

#define	PMPR_OUT	plhs[0]





/*  macros */

#define	max(A, B)	((A) > (B) ? (A) : (B))

#define	min(A, B)	((A) < (B) ? (A) : (B))





/*  function prototypes */



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



/*  Internal functions */





/*

---------------------------------------------------------------------------- 

    Compute the offset which points to a first/second order term in the 

    PmPr matrix.  i and j are the lags in the second-order term, k is

    the lag in the first-order term.



    Calculating the offset required for column (ij) is complicated by 

    the structure of the regressor matrix, which does not include any

    redundant columns.  Thus, the order of second order terms is

    (0,0) (0,1) ...(0,lags-1), (1,1) (1,2), ... (1,lags-1), (2,2) ect...

    notice that (1,0) is not present in the above... 

----------------------------------------------------------------------------

*/ 



#ifdef __STDC__

int offset21(

             int             i,

             int             j,

             int             k,

             int             lags,

             int             numpar

	)

#else

int offset21(i,j,k,lags,numpar)

             int                 i;

             int                 j;

             int                 k;

             int                 lags;

             int                 numpar;

#endif

{

  int numcols, offset;





  offset = (numpar+1+lags) + k*numpar;

  numcols = j-i + i*lags - i*(i-1)/2;

  offset += numcols;

  return offset;

}  





/*

---------------------------------------------------------------------------- 

    Compute the offset which points to a second/second order term in the 

    PmPr matrix.  i and j are the lags in the first second-order term, 

    k and l are the lags in the second second-order term.

----------------------------------------------------------------------------

*/ 



#ifdef __STDC__

int offset22(

             int             i,

             int             j,

             int             k,

             int             l,

             int             lags,

             int             numpar

	)

#else

int offset22(i,j,k,l,lags,numpar)

             int                 i;

             int                 j;

             int                 k;

             int                 l;

             int                 lags;

             int                 numpar;

#endif

{

  int numcols, numrows,offset;





  offset = (numpar+lags+1) + lags*numpar;

  numcols = j-i + i*lags - i*(i-1)/2;

  numrows = l-k + k*lags - k*(k-1)/2;

  offset += (numpar*numrows) + numcols;

  return offset;

}  







/*

---------------------------------------------------------------------------- 

    Place the value proj in the ijk element of the matrix pointed to by

    pmpr.  Also place this value in the "symmetric" locations, if they

    exist.  Assumes that k <= i <= j

----------------------------------------------------------------------------

*/ 









#ifdef __STDC__

void put21(

             double          pmpr[],

             double          proj,

             int             i,

             int             j,

             int             k,

             int             lags,

             int             numpar

	)

#else

void put21(pmpr,proj,i,j,k,lags,numpar)

             double              pmpr[];

             double              proj;

             int                 i;

             int                 j;

             int                 k;

             int                 lags;

             int                 numpar;

#endif

{

  *(pmpr + offset21(i,j,k,lags,numpar)) = proj;

  if (i > k)

    *(pmpr + offset21(k,j,i,lags,numpar)) = proj;

  if (j > i)

    *(pmpr + offset21(k,i,j,lags,numpar)) = proj;

}  







/*

---------------------------------------------------------------------------- 

    Place the value proj in the ijkl element of the matrix pointed to by

    pmpr.  Also place this value in the "symmetric" locations, if they

    exist.  Assumes k <= l <= i <= j

----------------------------------------------------------------------------

*/ 



#ifdef __STDC__

void put22(

             double          pmpr[],

             double          proj,

             int             i,

             int             j,

             int             k,

             int             l,

             int             lags,

             int             numpar

	)

#else

void put22(pmpr,proj,i,j,k,l,lags,numpar)

             double              pmpr[];

             double              proj;

             int                 i;

             int                 j;

             int                 k;

             int                 l;

             int                 lags;

             int                 numpar;

#endif

{

  if(k>lags)

    printf("spot 1 %d %d %d %d \n", i,j,k,l);

  *(pmpr + offset22(i,j,k,l,lags,numpar)) = proj;

  if (i > l)

    *(pmpr + offset22(l,j,k,i,lags,numpar)) = proj;

  if (j > l)

    *(pmpr + offset22(l,i,k,j,lags,numpar)) = proj;  

}





/*

---------------------------------------------------------------------------- 

  take the upper triangular entries in the numpar by numpar matrix

  starting at pmpr, put symmetric elements under the diagonal.

----------------------------------------------------------------------------

*/ 



#ifdef __STDC__

void symmetrize(

         double              pmpr[],

             int             numpar

	)

#else

void symmetrize(pmpr,numpar)

             double              pmpr[];

             int                 numpar;

#endif



{

  int i,j;

  for(i=0;i<numpar;i++){

    for(j=i;j<numpar;j++){

      *(pmpr + numpar*j + i) = *(pmpr + numpar*i + j);

    }

  }

}







#ifdef __STDC__

void fill_pmpr(

         double              pmpr[],

             int             numpar

	)

#else

void fill_pmpr(pmpr,numpar)

             double              pmpr[];

             int                 numpar;

#endif



{

  int i,j;

  for(i=0;i<numpar;i++){

    for(j=i;j<numpar;j++){

      *(pmpr + numpar*j + i) = 0;

      *(pmpr + numpar*i + j) = 0;

    }

  }

}





/*

----------------------------------------------------------------------------

  Function which calculates the first row of the PmPr matrix, which

  contains projections of regressors on a column of 1s. i.e. the zero

  order projections.

----------------------------------------------------------------------------

*/



#ifdef __STDC__

void order_0(

             double              pmpr[],

	     double	         phi1[],

	     double              u[],

	     int                 lags, 

	     unsigned int        data_length

	)

#else

void order_0(pmpr,phi1,u,lags,data_length)

             double              pmpr[];

	     double	         phi1[];

	     double	         u[];

             int                 lags;

	     unsigned int        data_length;

#endif

{

  double mu,correction,proj,*ptr_i,*ptr_j;

  int    i,j,offset;





/*  calculate input mean */

  mu = 0;

  for (i=0;i<data_length;i++) { 

    mu += *(u+i);

  }

  mu = mu/data_length;



/*  

    Calculate first column, projection of regressors onto

    a constant vector --> first element is simply 1 (projection

    of a constant vector onto itself)            

*/



  *pmpr = 1;



/*

  next, we have projections of lagged inputs onto the constant vector

  therefore, we get the input mean (0 lag), then the mean minus those

  points which "fall off the end" of the projection, when the input

  gets lagged progressively more and more.

*/



  proj = mu;

  ptr_i = u + data_length-1;               /* last element in u */

  *(pmpr+1) = proj;                        /* first term is input mean */ 

  for (i=1;i<lags;i++){

    proj -= *(ptr_i - i + 1)/data_length;  /* remove last remaining */

    *(pmpr+1+i) = proj;                    /* input value from term */

    }





/* 

  here, we have the projection of second-order regressors,  i.e. pairs

  of lagged inputs, onto the constant vector.  Corrections are a

  little more complicated, because adjacent terms are no longer shifted

  versions of each other, due to the ordering of the second-order

  regressors.



  The first 2/0 element is the auto-correlation a 0 lag.  The next 

  (lags-2) elements are the auto-correlation at lags of 1 through 

  (lags-1).  The next element contains the 0 lag auto-correlation,

  but delayed one sample, so that u(N)u(N)/N falls off the end.  

  Similarly, the next (lags-3) terms contain the auto-correlation 

  at lags of 1 through (lags-2), with one correction.



  We start with an uncorrected element (one of the frist lags entries),

  after 1 correction, it appears "lags" elements later, after the 

  second correction, it appears (lags-1) elements farther on.  When

  the size of the jump is smaller than the delay in the autocorrelation

  element, we have reached the end, and must start with a new 

  autocorrelation element.  

  

*/







  for(i=0;i<lags;i++){

    offset = lags + 1 + i;

    proj = *(phi1+i);

    *(pmpr + offset) = proj;

    ptr_i = u + data_length - 1;

    ptr_j = ptr_i - i;

    for(j=lags;j>i+1;j--){

      offset += j;

      correction = *ptr_i-- * *ptr_j--;

      proj -= correction/data_length;

      *(pmpr + offset) = proj;

    }

  }

}



/*

----------------------------------------------------------------------------

  Function which calculates the second diagonal block of the PmPr matrix, 

  which contains projections of first-order regressors onto other first-order 

  regressors.  These are all based on corrected versions of the input 

  autocorellation function.

----------------------------------------------------------------------------

*/



#ifdef __STDC__

void order_11(

             double              pmpr[],

	     double	         phi1[],

	     double              u[],

	     int                 lags,

             int                 numpar, 

	     unsigned int        data_length

	)

#else

void order_11(pmpr,phi1,u,lags,numpar,data_length)

             double              pmpr[];

	     double	         phi1[];

	     double	         u[];

             int                 lags;

             int                 numpar;

	     unsigned int        data_length;

#endif

{

  int      i, j, offset, new_offset;

  double   proj, correction, *ptr_i, *ptr_j;

  

/*

  Entries corresponding to pairs of first-order terms.  The idea is to

  move diagonally down the PmPr matrix, since all terms on a given diagonal 

  are based on the same cross-correlation entry, since the relative lags are 

  the same (in fact, at each step, only one additional product falls off the 

  end of the correlation and needs correction.

*/



  for(i=0;i<lags;i++) {

    offset = numpar + i + 1;  	     /* offset by numpars moves down one */

    *(pmpr+offset) = *(phi1+i);     /* row, + 1 to get to diagonal */ 

    ptr_i = u + data_length - i -1;

    ptr_j = u + data_length - 1;  	

    for(j=1;j<lags-i;j++) {

      correction = *ptr_i-- * *ptr_j-- / data_length;

      new_offset = offset + numpar + 1;  		

      *(pmpr + new_offset) = *(pmpr+offset) - correction;

      offset = new_offset;  		

    }

  }

}





/*

----------------------------------------------------------------------------

  Function which calculates the second/third off-diagonal block of the 

  PmPr matrix,  which contains projections of first-order regressors onto 

  second-order regressors.  These are all based on corrected versions of the 

  second (somethimes called third) order input autocorellation function.

----------------------------------------------------------------------------

*/



#ifdef __STDC__

void order_12(

             double              pmpr[],

	     double	         phi2[],

	     double              u[],

	     int                 lags,

             int                 numpar, 

	     unsigned int        data_length

	)

#else

void order_12(pmpr,phi2,u,lags,numpar,data_length)

             double              pmpr[];

	     double	         phi2[];

	     double	         u[];

             int                 lags;

             int                 numpar;

	     unsigned int        data_length;

#endif

{

  int i,j,k;

  double proj,correction,*ptr_k,*ptr_j,*ptr_i;



  for(i=0;i<lags;i++){

    for(j=i;j<lags;j++){

      k = 0;

      proj = *(phi2 + lags*i + j);

      put21(pmpr,proj,i,j,k,lags,numpar);

      ptr_k = u + data_length - 1;

      ptr_i = ptr_k - i;

      ptr_j = ptr_k - j;

      for(k=1;k<(lags-j);k++){

        correction = (*ptr_k--) * (*ptr_j--) * (*ptr_i--);

        proj -= correction/data_length;

        put21(pmpr,proj,i+k,j+k,k,lags,numpar);

      }

    }

  }

}





/*

----------------------------------------------------------------------------

  Function which calculates the third diagonal block of the PmPr matrix,  

  which contains projections of second-order regressors onto other 

  second-order regressors.  These are all based on corrected versions of 

  the third (somethimes called fourth) order input autocorellation function.

----------------------------------------------------------------------------

*/



#ifdef __STDC__

void order_22(

             double              pmpr[],

	     double	         phi3[],

	     double              u[],

	     int                 lags,

             int                 numpar, 

	     unsigned int        data_length

	)

#else

void order_22(pmpr,phi3,u,lags,numpar,data_length)

             double              pmpr[];

	     double	         phi3[];

	     double	         u[];

             int                 lags;

             int                 numpar;

	     unsigned int        data_length;

#endif

{

  int i,j,k,l;

  double  proj,correction, *ptr_i, *ptr_j, *ptr_k, *ptr_l;

  

  for(l=0;l<lags;l++){

    for(i=l;i<lags;i++){

      for(j=i;j<lags;j++){

        k=0;

        proj = *(phi3 + matrix_ref(lags,i,j,l));

        put22(pmpr,proj,i,j,k,l,lags,numpar);

        ptr_k = u + data_length - 1;

        ptr_i = ptr_k - i;

        ptr_j = ptr_k - j;

        ptr_l = ptr_k - l;

        for(k=1;k<lags-j;k++){

          correction = (*ptr_i--)*(*ptr_j--)*(*ptr_k--)*(*ptr_l--);

          proj -= correction/data_length;

          put22(pmpr,proj,i+k,j+k,k,l+k,lags,numpar);

	}

      }

    }

  }  

}





/*

    This is the main computation routine. 



*/



#ifdef __STDC__

void fast_pmpr(

         double              pmpr[],

	     double	         phi1[],

	     double	         phi2[],

	     double          phi3[],

	     double          u[],

	     int             lags, 

             int             numpar,

	     unsigned int    data_length

	)

#else

void fast_pmpr(pmpr,phi1,phi2,phi3,u,lags,numpar,data_length)

             double              pmpr[];

	     double	         phi1[];

	     double	         phi2[];

	     double	         phi3[];

	     double	         u[];

             int                 lags;

	     int                 numpar; 

	     unsigned int        data_length;

#endif





              

{

  fill_pmpr(pmpr,numpar);

  order_0(pmpr,phi1,u,lags,data_length);

  order_11(pmpr,phi1,u,lags,numpar,data_length);

  order_12(pmpr,phi2,u,lags,numpar,data_length);

  order_22(pmpr,phi3,u,lags,numpar,data_length);

  symmetrize(pmpr,numpar);

}



/*   end of function fast_pmpr    */





/*

    This is the entry point for MATLAB



*/



#ifdef __STDC__

void mexFunction(

		 int		nlhs,

		 mxArray        *plhs[],

		 int		nrhs,

		 const mxArray	*prhs[]

		 )

#else

void mexFunction(nlhs, plhs, nrhs, prhs)

int nlhs, nrhs;

mxArray *plhs[];

const mxArray *prhs[];

#endif



{

  double           *phiu2,*phiu3,*phiu4;

  double           *u, *pmpr;

  unsigned int     lags, data_length, numpar,m, n;

 



  /* Check for proper number of arguments */

  

  if (nrhs == 0) {     /*  print help message */

    printf("\n MEX-file implementation of fast_pmpr.m \n");

    printf("\n Usage:  pmpr = fast_pmpr(phiu2,phiu3,phiu4,u); \n");

    printf("\t phiu2:  first-order cross-correlation \n");

    printf("\t phiu3:  second-order cross-correlation \n");

    printf("\t phiu4:  third-order cross-correlation \n");

    printf("\t u:  input signal \n");

  }

  else if (nrhs != 4) {

    mexErrMsgTxt("fast_pmpr requires four input arguments.");

  } else if (nlhs > 1) {

    mexErrMsgTxt("fast_pmpr requires one output argument.");

  } else {   /* proceed with operator */    



  /* Assign pointers to the various parameters */



    phiu2 = mxGetPr(P1_IN);

    phiu3 = mxGetPr(P2_IN);

    phiu4 = mxGetPr(P3_IN);

    u = mxGetPr(U_IN);

  

    m = mxGetM(U_IN);

    n = mxGetN(U_IN);

    data_length = max(m,n);

    m = mxGetM(P1_IN);

    n = mxGetN(P1_IN);

    lags = max(m,n);



  /* Create a matrix for the return argument */



    numpar = 1 + lags + lags*(lags+1)/2;

    PMPR_OUT = mxCreateDoubleMatrix(numpar, numpar, mxREAL);

    pmpr = mxGetPr(PMPR_OUT);





  /* Do the actual computations in  other subroutines */



    fast_pmpr (pmpr,phiu2,phiu3,phiu4,u,lags,numpar,data_length);

    return;

  }

}





















