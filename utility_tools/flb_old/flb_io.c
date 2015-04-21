/* flb_io.c 

26 aug 98 REK Changes for V5.2 and MS C

*/


#  include "stdlib.h"
#  include "stdio.h"
#  include "string.h"
#  include "filelb.h"
#  include "debug.h"

	void read_vms_flb_details (FILE *flb_ptr, FLB_DETAILS *details);

void write_flb_details (flb_ptr, details)
FILE *flb_ptr;
FLB_DETAILS *details;
{

  QUARTET string_len, vnum;
  INT_T i, num_names;

  vnum = 2;
  fwrite (&vnum, sizeof(QUARTET), 1, flb_ptr);
  fwrite (&(details->nchan), sizeof(QUARTET), 1, flb_ptr);
  fwrite (&(details->length), sizeof(QUARTET), 1, flb_ptr);
  fwrite (&(details->format), sizeof(QUARTET), 1, flb_ptr);
  string_len = strlen (details->domain_name);
  fwrite (&string_len, sizeof(QUARTET), 1, flb_ptr);
  fwrite ((details->domain_name), sizeof(char), string_len, flb_ptr);
  fwrite (&(details->increment), sizeof(float), 1, flb_ptr); 
  fwrite (&(details->domain_start), sizeof(float), 1, flb_ptr);
  string_len = strlen (details->comment);
  fwrite (&string_len, sizeof(QUARTET), 1, flb_ptr);
  fwrite ((details->comment), sizeof(char), string_len, flb_ptr);
  if (details->nchan > FLB_CHAN_MAX) {
    num_names = 1;
    }
  else {
    num_names = details->nchan;
  }
  for (i=0; i < num_names; i++){
    string_len = strlen (details->name[i]);
    fwrite (&string_len, sizeof(QUARTET), 1, flb_ptr);
    fwrite ((details->name[i]), sizeof(char),string_len, flb_ptr);
  }
  fwrite ((details->min), sizeof (double), num_names, flb_ptr);
  fwrite ((details->max), sizeof (double), num_names, flb_ptr);
  return;
}



INT_T read_flb_details ( flb_ptr, details )


     FILE  *flb_ptr;
     FLB_DETAILS *details;
{

  INT_T err_flag = -1;
  QUARTET string_len;
  INT_T i, num_names;
  if(fread (&(details->version_num), sizeof(QUARTET), 1, flb_ptr) != 1) {
  DEBUG (
	DEBUG_LOC("read_flb_details");
	 printf ("eof\n");
 )
    return err_flag;
     };

  DEBUG(
	printf("\n");
  	DEBUG_LOC("read_flb_details");
        printf("details->version_num = %d\n\n" COMMA details->version_num); 
  )

  if (details->version_num == 1) {
    read_vms_flb_details ( flb_ptr, details );
  }
  else if (details->version_num ==2){
    fread (&(details->nchan), sizeof(QUARTET), 1, flb_ptr);
    fread (&(details->length), sizeof(QUARTET), 1, flb_ptr);
    fread (&(details->format), sizeof(QUARTET), 1, flb_ptr);
    fread (&string_len, sizeof(QUARTET), 1, flb_ptr);
    if (string_len > FLB_NAME_MAX) {
      printf ("Domain name too long: %s\n", details->domain_name);
      return -1;
    }

    fread ((details->domain_name), sizeof(char), string_len, flb_ptr);
    details->domain_name[string_len] = '\0';
    fread (&(details->increment), sizeof(float), 1, flb_ptr); 
    fread (&(details->domain_start), sizeof(float), 1, flb_ptr);
    fread (&string_len, sizeof(QUARTET), 1, flb_ptr);
    if (string_len > FLB_COMMENT_MAX) {
      printf ("Comment too long: %s\n", details->comment);
      return -1;
    }
    fread ((details->comment), sizeof(char), string_len, flb_ptr);
    details->comment[string_len]= '\0';
    if (details->nchan > FLB_CHAN_MAX) {
      num_names = 1;
      }
    else {
      num_names = details->nchan;
	}
    for (i=0; i < num_names; i++){
      fread (&string_len, sizeof(QUARTET), 1, flb_ptr);
      if (string_len > FLB_NAME_MAX) {
	printf ("Channel name too long: %s\n", details->name[i]);
        return -1;
      }
      fread ((details->name[i]), sizeof(char), string_len, flb_ptr);
    details->name[i][string_len] = '\0';
    }
    fread ((details->min), sizeof (double), num_names, flb_ptr);
    fread ((details->max), sizeof (double), num_names, flb_ptr);
  }
  else {
    printf ("Bad FLB version number\n");
    return -1;
  }
  return 1;
}




void read_vms_flb_details ( flb_ptr, details )

     FILE  *flb_ptr;
     FLB_DETAILS *details;
{
  INT_T i, count;

  fread (&(details->nchan), sizeof(QUARTET), 1, flb_ptr);
  fread (&(details->length), sizeof(QUARTET), 1, flb_ptr);
  if (details->nchan > 1 ) {
    fseek ( flb_ptr, 4*(details->nchan -1), SEEK_CUR);
  }
  fread (&(details->format), sizeof(QUARTET), 1, flb_ptr);
  if (details->nchan > 1 ) {
    fseek ( flb_ptr, 4*(details->nchan -1), SEEK_CUR);
  }
  fread ((details->domain_name), 8, 1, flb_ptr);
  details->domain_name[8]= '\0';
  if (details->nchan > 1 ) {
    fseek ( flb_ptr, 8*(details->nchan -1), SEEK_CUR);
  }
  fread (&(details->increment), sizeof(float), 1, flb_ptr); 
  if (details->nchan > 1 ) {
    fseek ( flb_ptr, 4*(details->nchan -1), SEEK_CUR);
  }

  fread (&(details->domain_start), sizeof(float), 1, flb_ptr);
   if (details->nchan > 1 ) {
    fseek ( flb_ptr, 4*(details->nchan -1), SEEK_CUR);
  }

  fseek ( flb_ptr, 4* details->nchan, SEEK_CUR);
  for (i=0; i< details->nchan; i++){
    fread (details->name[i], 1,16, flb_ptr);
    details->name[i][16] = '\0';
  }
  fread (details->comment, sizeof(char), 80, flb_ptr);
  details->comment[79]='\0';
  count = 88 + (details->nchan)*44;
  fseek ( flb_ptr, 512-count, SEEK_CUR);
  return;
}
















