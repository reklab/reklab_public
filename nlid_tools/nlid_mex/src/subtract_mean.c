/*

  subtract_mean.c
  This function takes a pointer to an input vector, x, 
  and its length, data_length 
  The mean of x is calculated and subtracted from x.		 

  subtract_mean (*x,data_length)

  David Westwick   August 23, 1993
 
*/


#ifdef __STDC__
void subtract_mean (
      double *x,				  /*  input channel  */
      unsigned int  data_length			 /* length of chanel */
    )
#else
void subtract_mean(x,data_length)
      double *x;
      unsigned int data_length;
#endif
{
  int i;
  double sum, *dummy_ptr;
  
  sum = 0.0;
  dummy_ptr = x;
  for (i=0;i<data_length;i++) {
    sum += *dummy_ptr++;
    }
  dummy_ptr = x;
  sum = sum/data_length;
  for (i=0;i<data_length;i++) {
    *dummy_ptr++ -= sum;
    } 	 
  }  /* end of function subtract_mean */

