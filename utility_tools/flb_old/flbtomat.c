/* 
     flbtomat - mex file convert flb files to mat files
     
Options:
[  x, comment, names, domain, incr, start] = flbtomat (file, case, option)

13 May 92  - close flb_ptr before calling mex_error.
14 May 92  - corrected close error - use fclose rather than close
18 Dec 92  - modify output format
25 May 93  - add bitmap -to- binary matrix extraction option
             (see mattoflb.c for complementary half) [jd.]
4  Aug 93  - Machine-independent version, and porting to MSDOS [jd.]
             (addition of "flbport.h", which defines machine-indep word sizes)
6 Jan 94   - modified to use standard calls for matlab 4.1 for use with alpha
15 May 97  - REK modify for Matlab5 
18 Nob 97  - Problem with increment
26 Aug 98  - V05-03 Many changes to make standard.
10 Sep 98  - V05-04 Fixed memory problem with use of realloc. Need to realloc
			 with the current length (pcurrent) + the length of the 
			 comment (details.comment) + 7 ("%3i: " *AND* the "|").
*/

//#ifdef MSDOS			/* different headers for MSDOS */
#  include <stdio.h>
#  include <mex.h>
#  include <math.h>
#  include <string.h>
#  include "filelb.h"
#  include "debug.h"

   /* forward declarations, for TurboC... */

   extern INT_T read_vms_flb_details ( FILE *, FLB_DETAILS * );
   extern INT_T read_flb_details ( FILE *, FLB_DETAILS *); 
   extern void write_flb_details ( FILE *, FLB_DETAILS *); 
   extern loadmat( FILE *,INT_T *,char *,INT_T *,INT_T *,INT_T *,double **,double **);
   extern savemat( FILE *,INT_T,char *,INT_T,INT_T,INT_T,double *,double *);
   extern INT_T load_header(FILE *,INT_T *,char *,INT_T *,INT_T *,INT_T *,INT_T *);

   INT_T read_flb_case (FLB_DETAILS *, FILE *);
   INT_T read_vms_flb_case (FLB_DETAILS *, FILE *);
   void skip_flb_case (FLB_DETAILS *, FILE *);
   void skip_vms_flb_case (FLB_DETAILS *, FILE *);
   void print_flb_header (FLB_DETAILS *, INT_T, INT_T);
   void flbtomat_help (void); 


/*
#else               /* headers for Unix, VMS, etc...

#  include "stdio.h"
#  include <stdlib.h>
#  include "math.h"

#  include "string.h"
#  include "debug.h"
#  include "mex.h"
#  include "filelb.h"

   extern INT_T read_flb_details ( FILE *, FLB_DETAILS *); 

#endif
*/

#define FLB_NAME prhs [0]
#define CASE_NUM prhs [1]
#define OPTION   prhs [1]
#define Y_OUT plhs [0]
#define Y_COMMENT plhs [1]
#define Y_INCREMENT plhs [2]
#define Y_DOMAIN    plhs [3]
#define Y_NAMES     plhs [4]
#define Y_START     plhs [5]

/* Some defines used for bitmap conversion */
#define BITMAP	(-1)	/* for details.format identifier */
			/* == -sizeof(BYTE) by the way!  */
#define BYTE	unsigned char
#define LSBit	1
#define	MSBit	(1<<7)

unsigned QUARTET case_num; 
char vname[80];
float *flb_data;   /*  Array to hold data points */
char   *flb_data_char;
short  *flb_data_i2;
unsigned short int *flb_data_unsigned;
BYTE   *flb_data_byte;    
double *mat_data;

/* ------------------------------------------------------------------------
   flbtomat - main routine
   ------------------------------------------------------------------------
*/

void mexFunction (nlhs, plhs, nrhs, prhs)
     INT_T  nlhs;
     INT_T nrhs;
     mxArray *plhs[];
     const mxArray *prhs[];

{
  FILE *flb_ptr;
  FLB_DETAILS details;
  char *flb_file;
  char *pcurrent;
  char flb_file_name[80], option [16], case_char[16];
  double scale, offset;
  unsigned QUARTET input_case_num;
  INT_T i,itemp, length;
  INT_T  type, mrows, ncols;
  INT_T col,row,flb_index, mat_index; 
  INT_T flb_flag, index_flag, display_flag, length_flag, output_flag;
  INT_T  verbose_flag;
  int n;
  size_t lnew;
  double *increment, *p;
  input_case_num = MAX_VAL(unsigned QUARTET);
  display_flag = 0;
  index_flag = 0;
  verbose_flag = 2;
  output_flag = 1;
  length_flag =0; 

 if (nrhs == 0) {
    flbtomat_help ();
    return;
  }

  if (nlhs > 6 ) {
    mexErrMsgTxt ("Too many output channels specified\n");
    }

/*  Open filelb channel file */  	

  flb_file = flb_file_name;
  n=mxGetN(FLB_NAME)+1;
  mxGetString (FLB_NAME,flb_file,  n);
    if((flb_ptr=fopen (flb_file, "rb"))== NULL){
      mexErrMsgTxt("Error opening FLB file");
   }

  /*  Case number */
    
  if (nrhs == 1) {
    display_flag = 1;
    output_flag = 0; 
    input_case_num = MAX_VAL(unsigned QUARTET);
    length_flag = 1;
    index_flag =1;
  }
  
  if (nrhs == 2) {

      if ( mxIsChar(OPTION)) {
	index_flag = 1;
	n=mxGetN(OPTION)+1;
	mxGetString(OPTION, option, 2);
	if ( strcmp(option, "v")==0) {
	  verbose_flag = 1;
	  display_flag = 1;
	  output_flag = 0;
	  length_flag=1;
	  input_case_num = MAX_VAL(unsigned QUARTET);		
	}
	else if ( strcmp(option,"q") ==0) {
	  verbose_flag = 0;
	  display_flag = 0;
	  output_flag =0;
	  length_flag=1;
	  input_case_num = MAX_VAL(unsigned QUARTET);		
	}
	else {
	mexErrMsgTxt("flbtomat: bad option value");
      }
      }
      else {
	input_case_num = ( INT_T ) mxGetScalar(CASE_NUM);
	index_flag =0;
      }
  }

    DEBUG(   DEBUG_LOC("user_fcn");
	)

  /* Read case header */

  DEBUG(   DEBUG_LOC("user_fcn");
	   mexPrintf("input_case_num = %u\n", input_case_num);
	)

  case_num = 0;
/* Initialize pointers for output of comment strings */
  if ((nlhs>1) & (index_flag==1)){
    lnew=0;
  }
  else {
    index_flag = 0;
  }

  while ((case_num<input_case_num)){
    case_num++;

  DEBUG(
	mexPrintf("\n");
	DEBUG_LOC("user_fcn");
	mexPrintf("case_num = %d\n", case_num);
	mexPrintf("Read details\n");
  )

     flb_flag=read_flb_details (flb_ptr, &details);
    if (flb_flag ==-1){
      if ((output_flag == 1) & (length_flag == 0)){
	fclose(flb_ptr);
	mexErrMsgTxt ("Error reading flb file \n");
      }
      if (length_flag == 1){
	Y_OUT=mxCreateDoubleMatrix (1,1,mxREAL);
	p = mxGetPr (Y_OUT);
	*p= (double) (case_num-1);
      }
      if (index_flag == 1){
	Y_COMMENT = mxCreateString (pcurrent);
	mxFree(pcurrent);
      }
     else {
      }
      
      return;
    }

    if (details.nchan > FLB_CHAN_MAX) {
           mexPrintf("flbtomat: FLB file has many channels. Treating as an enesemble\n");
	   }   

    if (display_flag == 1) {
      print_flb_header(&details, case_num, verbose_flag);
      if (verbose_flag == 2) {
	verbose_flag =3;
      }
    }

    if (index_flag == 1) {
      sprintf(case_char, "%3i: ",case_num);
      if (lnew ==0){
	lnew = strlen (details.comment) + 6;
	pcurrent = mxCalloc (lnew,1);
	strcpy (pcurrent,case_char);
	strcat (pcurrent, details.comment);
      }
      else {
	lnew = strlen (pcurrent) + strlen (details.comment) + 7;
	pcurrent = mxRealloc (pcurrent,lnew);
	strcat ( pcurrent,"|");
	strcat (pcurrent,case_char);
	strcat (pcurrent, details.comment);
	
      }
   
    }
   
    
    /* Read data and convert */
    if ((output_flag ==1 ) & (case_num == input_case_num)) {
      if (read_flb_case(&details, flb_ptr)==0){
	fclose (flb_ptr);
	mexErrMsgTxt ("Error reading data\n");
      }
    /* Convert to MATLAB double precision */
      length = details.length;
      mrows = length;
      ncols = details.nchan;
      Y_OUT = mxCreateDoubleMatrix ( mrows, ncols, mxREAL);
      mat_data=mxGetPr(Y_OUT);
      flb_index = 0;		/* index to step thru flb buffers */

      if (details.format == BITMAP){  
      /* Treat bitmapped data separately... */
      /* (because number of entries == number of BITS, not bytes) */

        /* first declare & define some very local variables */
        unsigned long  size = mrows*ncols,     /* know when to stop */
	i = 0;               /* BIT index */
	
        BYTE bit,         		/* (mask) to check if 1 or 0 */
	byteme;                    /* buffer to hold current byte */

      /* this loop examines a byte at a time, and within that byte */
      /* looks at each bit in turn (by masking out the others)...  */
      /* Then decides whether the matrix entry should be a 1 or 0  */
      /* depending on whether the bit is set/reset.                */

        while( i<size ){
          byteme = flb_data_byte[flb_index++]; 	/* get new byte */
          bit = (BYTE)LSBit;			/* mask out first bit */
            
          while( bit && (i<size) ){   /* still on this bit and have data */
            mat_data[i++] = (double)( (byteme & bit) ? 1:0 ); 
            bit <<= 1; 			/* mask out next bit */
          }
        }
      }  /* end of BITMAP case */

      else{
      /* all other types of data are stored in multiples of one byte */

        for (col=0;col<ncols;col++) {
	  scale=(details.max[col]-details.min[col])/UNSIGNED_SHORT_MAX;
	  offset = details.min[col];
	  for (row=0;row<mrows;row++) {
	    mat_index = col*mrows + row;
	    if (details.format == 1) {
	      itemp  = *(flb_data_char + flb_index);
	      *(mat_data+mat_index) = (double) itemp;
	    }    
	    else if (details.format == 2) {
	      itemp  = *(flb_data_i2 + flb_index);
	      *(mat_data+mat_index) = (double) itemp;
	    }
	    else if (details.format == -2) {
	      itemp  = *(flb_data_unsigned + flb_index);
	      *(mat_data+mat_index) = ((double) itemp)*scale + offset;
	    }  
	    else if (details.format == 4) {
	      *(mat_data+mat_index) = *(flb_data + flb_index);
	    }             
	  flb_index++;;
	  } /* for row */
        }  /* for col */
      } /* end of else (not BITMAP) case */
      
      type =0;
      if (details.format ==1) mxFree (flb_data_char);
      if (details.format ==2) mxFree (flb_data_i2);
      if (details.format ==-2) mxFree (flb_data_unsigned);
      else if (details.format ==4) mxFree (flb_data);
    } /* if output_flag==1 and case_num==input_case_num */

    else {
      skip_flb_case(&details, flb_ptr) ;
    }
  }

  /* Comment */
   if (nlhs >= 2) {
       Y_COMMENT = mxCreateString (details.comment);
   }
/*  Domain increment */

  if (nlhs >= 3) {
    Y_INCREMENT = mxCreateDoubleMatrix(1,1, mxREAL);
    increment = mxGetPr(Y_INCREMENT);
    *increment = (double) details.increment;
  }
  /* Domain name */
  if (nlhs >= 4) {
    Y_DOMAIN = mxCreateString (details.domain_name);
    }
  /* Variable names */
  if (nlhs >= 5) {
    pcurrent = mxCalloc ( FLB_NAME_MAX * FLB_CHAN_MAX, sizeof (char));
    for (i=0; i < details.nchan; i++) {
      if (i==0) {
	strcpy (pcurrent, details.name[i]);
      }
      else {
	strcat (pcurrent,"|");
	strcat (pcurrent,details.name[i]);
      }
    }
    Y_NAMES = mxCreateString(pcurrent);
	mxFree(pcurrent);
  }
  /* Domain start */

  if (nlhs >= 6) {
    Y_START = mxCreateDoubleMatrix( 1,1,  mxREAL);
    p = mxGetPr (Y_START);
    *p=details.domain_start;
  }

/* End of file - close files and exit */
  fclose(flb_ptr);
  
}


/* ---------------------------------------------------------------------
   read_flb_case
  -----------------------------------------------------------------------*/

INT_T read_flb_case (details, flb_ptr)
     FLB_DETAILS *details;
     FILE *flb_ptr;

{

  INT_T nitems, nbytes, cf, cf_size;
  INT_T nsamp;
  if (details->version_num == 1 ) {
    nitems = read_vms_flb_case (details, flb_ptr);
  }
  else {
    cf = details->format;
    cf_size = cf;
    if (cf_size <0) cf_size = - cf_size;
    nsamp = details->nchan * details->length ;
    nbytes = nsamp *cf_size;
    if ( cf == 1) {
      flb_data_char = (char *) mxCalloc ( nsamp, cf_size );
      nitems=fread(flb_data_char, cf_size, nsamp, flb_ptr);
    }
    else if ( cf==2) {
      flb_data_i2 = (short *) mxCalloc ( nsamp, cf_size );
      nitems= fread(flb_data_i2, cf_size, nsamp, flb_ptr);
    }
    else if (cf==-2) {
      flb_data_unsigned=(unsigned short *) mxCalloc (nsamp, cf_size );
      nitems= fread(flb_data_unsigned, cf_size, nsamp, flb_ptr);
    }
    else if ( cf == 4) {
      flb_data = (float *) mxCalloc ( nsamp, cf_size );
      nitems= fread(flb_data, cf_size, nsamp,flb_ptr);
    }
    else if ( cf == BITMAP ){
       /* do a ceil(nsamp/8) because nsamp <--> number of BITS not bytes */
       nsamp = (INT_T)ceil((double)nsamp/8);

       flb_data_byte = (BYTE *) mxCalloc ( nsamp, cf_size );
       nitems=fread(flb_data_byte,cf_size, nsamp, flb_ptr);
    }
    
  }
  
  if (nitems != nsamp) {
    fclose(flb_ptr);
	mexPrintf ("Error reading input file\n");
  }
  return nitems;	
} /* read_flb_case () */


/* ---------------------------------------------------------------------
   read_vms_flb_case
  -----------------------------------------------------------------------*/

INT_T read_vms_flb_case (details, flb_ptr)
     FLB_DETAILS *details;
     FILE *flb_ptr;
{
  INT_T cf, cf_size, nsamp, nbytes, nitems, off, i, remainder;
  cf = details->format;
  cf_size = cf;
  if (cf_size <0) cf_size = - cf_size;
  nsamp =  details->length ;
  nbytes = nsamp *cf_size;
  remainder = 512 - nbytes%512;
  off = 0;
  if (cf_size == 2) {
    flb_data_i2 = (short *) mxCalloc ( nsamp* details->nchan, cf_size );
  }
  else if (cf_size == 4) {
    flb_data = (float *) mxCalloc( nsamp * details->nchan, cf_size);
  }
  for (i=0; i < details->nchan; i++){
    if ( cf == 1) {
      mexErrMsgTxt("Format 1 not supported for VMS FLB files\n");
    }
    else if ( cf==2) {
      nitems= fread(flb_data_i2+off, cf_size, nsamp, flb_ptr);
      if (remainder < 512) fseek (flb_ptr, remainder, SEEK_CUR);
    }
    else if (cf==-2) {
      mexErrMsgTxt("Format -2 not supported for VMS FLB files\n");
    }
    else if (cf==BITMAP){
      mexErrMsgTxt("Format BITMAP not supported for VMS FLB files\n");
    }
    else if ( cf == 4) {
      nitems= fread(flb_data+off, cf_size, nsamp,flb_ptr);
      if (remainder < 512) fseek (flb_ptr, remainder, SEEK_CUR);
    }
    off = off + nsamp;
  }
  return nitems;
}




/* ---------------------------------------------------------------------
   skip_flb_case
  -----------------------------------------------------------------------*/

void skip_flb_case (details, flb_ptr)
     FLB_DETAILS *details;
     FILE *flb_ptr;

{

  QUARTET cf;
  LONG_T nbytes, cf_size, nsamp;

  if (details->version_num == 1 ) {
    skip_vms_flb_case ( details, flb_ptr);
  }
  else {
    cf = details->format;
    cf_size = cf;
    if (cf_size <0) cf_size = - cf_size;
    nsamp = details->nchan * details->length ;
    nbytes = nsamp *cf_size;

    DEBUG(	DEBUG_LOC("skip_flb_case");
		mexPrintf("cf_size = %d\n", cf);
		mexPrintf("nsamp   = %d\n", nsamp);
		mexPrintf("nbytes  = %d\n", nbytes);
         )

    if ( cf == 1) {
      fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);
    }
    else if ( cf==2) {
           fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);
    }
    else if (cf==-2) {
           fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);    
    }
    else if (cf==BITMAP){
           nsamp = (INT_T)ceil((double)nsamp/8); /* divide byte count by 8 */
           fseek ( flb_ptr, cf_size*nsamp, SEEK_CUR);
    }
    else if ( cf == 4) {
      fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);
    }
  }
}
/* read_flb_case () */


/* ---------------------------------------------------------------------
   skip_vms_flb_case
  -----------------------------------------------------------------------*/

void skip_vms_flb_case (details, flb_ptr)
     FLB_DETAILS *details;
     FILE *flb_ptr;
{
  QUARTET cf;
  LONG_T cf_size, nsamp, nbytes, off, i, remainder;

  cf = details->format;
  cf_size = cf;
  if (cf_size <0) cf_size = - cf_size;
  nsamp =  details->length ;
  nbytes = nsamp *cf_size;
  remainder = 512 - nbytes%512;
  off = 0;
  for (i=0; i < details->nchan; i++){
    if ( cf == 1) {
      mexErrMsgTxt("Format 1 not supported for VMS FLB files\n");
    }
    else if ( cf==2) {
      fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);
      if (remainder < 512) fseek (flb_ptr, remainder, SEEK_CUR);
    }
    else if (cf==-2) {
      mexErrMsgTxt("Format -2 not supported for VMS FLB files\n");
    }
    else if (cf==BITMAP){
      mexErrMsgTxt("Format BITMAP not supported for VMS FLB files\n");
    }
    else if ( cf == 4) {
      fseek ( flb_ptr,  cf_size* nsamp, SEEK_CUR);
      if (remainder < 512) fseek (flb_ptr, remainder, SEEK_CUR);
    }
    off = off + nsamp;
  }
  return;
}

/* ---------------------------------------------------------------------- 
   print_flb_header
   ---------------------------------------------------------------------- */
void print_flb_header (details,case_num, flag )
     FLB_DETAILS *details;
     INT_T case_num, flag;
{
  INT_T i, num_cols; 
  if ( flag == 1 ) {
    mexPrintf("\nCase number = %d;",case_num);
    mexPrintf(" FLB Version: %d;",details->version_num);	
    mexPrintf(" Length: %d;",details->length);	
    mexPrintf(" Nchan: %d;",details->nchan);	
    mexPrintf(" Format: %d\n",details->format);	
    mexPrintf("Comment: %s \n",details->comment);
    mexPrintf("Domain Name: %s ;",details->domain_name);	
    mexPrintf(" Incr = %f ;",details->increment);
    mexPrintf(" Start= %f \n",details->domain_start);
    num_cols = details->nchan;
    if (num_cols > FLB_CHAN_MAX)
     num_cols =1;
    for (i=0; i< num_cols; i++) {
      mexPrintf("Channel %i : Name = %s ; Min = %f ; max = %f\n", i, 
	      details->name[i], details->min[i], details->max[i]);
    }
  }
  else if ( flag ==2) {
    mexPrintf("Case  NChan Length  Comment \n");
    mexPrintf("%4d %4d %6d;",case_num, details->nchan,details->length);
    mexPrintf(" %s \n",details->comment);
  }
  else if ( flag ==3) {
    mexPrintf("%4d %4d %6d;",case_num, details->nchan, details->length);
    mexPrintf(" %s \n",details->comment);
  }
}

void flbtomat_help (void) 
{

DEBUG_LOC("Unix Compilation Trial 2");	/* for use with "debug.h" */


  mexPrintf("\nflbtomat V05-04 REK\n");

mexPrintf("\n convert FILELB files to MATLAB format\n");
mexPrintf("\n  Usage:  [ x c i dn cn ds]= flbtomat(file, case/option)");
mexPrintf("\n\n outputs are: \n");
mexPrintf("\t x - data\n");
mexPrintf("\t c - comment\n");
mexPrintf("\t i - increment\n");
mexPrintf("\t dn - domain name\n");
mexPrintf("\t cn - channel name\n");
mexPrintf("\t ds - domain start\n");
mexPrintf("\n inputs are: \n");
mexPrintf("\t file - names of filelb file\n");
mexPrintf("\t case - case to recall\n");
mexPrintf("\t option - display option\n");
mexPrintf("\n Calling flbtomat with no parameters gives a help message\n");
mexPrintf("\n\n If flbtomat is called without specifying a case number,\n");
mexPrintf(" an index is displayed on the screen and the outputs are: \n");
mexPrintf("\n\t x - number of cases inthe file \n");
mexPrintf("\t c - concatenated comments from all cases \n");
mexPrintf("\t If option  = v a verbose index is displayed\n");
mexPrintf("\t If option  = q the index is not displayed\n");
}

