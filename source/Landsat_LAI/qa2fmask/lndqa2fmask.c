/**
 * !Descriptio
 * n
 * This program reads "pixel_qa" file in Landsat Collection 1 surface reflectance product and 
 * convert to compatible fmask file to the pre-collection Landsat fmask file
 *
 * !Input   Landsat surface reflectance pixel_qa file
 *          <bit num="0">fill</bit>
 *          <bit num="1">clear</bit>
 *          <bit num="2">water</bit>
 *          <bit num="3">cloud shadow</bit>
 *          <bit num="4">snow</bit>
 *          <bit num="5">cloud</bit>
 *          <bit num="6">cloud confidence</bit>
 *          <bit num="7">cloud confidence</bit>
 *          <bit num="8">unused</bit>
 *          <bit num="9">unused</bit>
 *          <bit num="10">unused</bit>
 *          <bit num="11">unused</bit>
 *          <bit num="12">unused</bit>
 *          <bit num="13">unused</bit>
 *          <bit num="14">unused</bit>
 *          <bit num="15">unused</bit>
 *          cloud confidence (bit 6,7): 00=none; 01-low; 10=medium; 11=high 
 *
 * !Output  unpacked QA bits and saves in binary format with value  
 *          0 => clear land pixel
 *          1 => clear water pixel
 *          2 => cloud shadow
 *          3 => snow
 *          4 => cloud
 *          255 => no observation
 *
 * !Developer: 
 *          Feng Gao (feng.gao@ars.usda.gov)
 *
 * !Revision:
 *
 * Original version - 06/2017
 * modified for Landsat Collection 1: 06/2017
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_STRLEN 1000
#define MAX_NCOLS  12000
#define FAILURE -1
#define SUCCESS 1

typedef struct {
  char fileName[MAX_STRLEN];  /* SR hdf file name */
  int nrows;           
  int ncols;
  float ulx;
  float uly;
} GRID_SR;    

void usage(char *);
int  readENVIHeader(GRID_SR *sensor);
void dec2bin(unsigned short int num, int *bitpattern, int limit);

/*#define DEBUG*/
#define DEBUG_irow 1421
#define DEBUG_icol 4046

int main(int argc, char *argv[])
{
  int i, irow, icol, acceptable_cc, bitp[16];
  unsigned char vmask, cc;
  unsigned short int c;
  char fname[MAX_STRLEN];

  GRID_SR lndsr;
  FILE *in, *out;

  /* should include 5/6 args: lndqa2fmask.exe [-lndsr] <lndsr_file> [-cmask] <out_mask> [acceptable_cc]*/ 
  if(argc != 5 && argc != 7) {
    usage(argv[0]);
    return FAILURE;
  }
  
  acceptable_cc = 0;
  /* parse command line */
  for(i=1; i<argc; i++){
    if(strcmp(argv[i],"-lndsr")==0)
      strcpy(lndsr.fileName, argv[++i]);
    else
      if(strcmp(argv[i],"-cmask")==0)
	strcpy(fname, argv[++i]);
      else
	if(strcmp(argv[i],"-acc")==0)
	  acceptable_cc = atoi(argv[++i]);
	else{
	  printf("\nWrong option:%s\n",argv[i]);
	  usage(argv[0]);
	  return FAILURE;
	}    
  } 
  
  if((in=fopen(lndsr.fileName, "rb"))==NULL) {
    printf("Open file %s error\n", lndsr.fileName);
    return FAILURE;
  }
    
  /* open file for output */
  if((out=fopen(fname, "wb"))==NULL) {
    printf("Open file %s error\n", fname);
    return FAILURE;
  }

  readENVIHeader(&lndsr);

  for(irow=0; irow<lndsr.nrows; irow++) 
    for(icol=0; icol<lndsr.ncols; icol++) {
      fread(&c, sizeof(short int), 1, in);  
      dec2bin(c, bitp, 16);    
      vmask = 255;
      cc = bitp[6]*2+bitp[5];
      if(bitp[1]==1) vmask = 0;
      else if(bitp[2]==1) vmask = 1;
      else if(bitp[3]==1) vmask = 2;
      else if(bitp[4]==1) vmask = 3;
      else if(bitp[5]==1) { 
	if(cc>=acceptable_cc) vmask = 4; 
	else vmask = 0;
      }
      fwrite(&vmask, 1, 1, out);      
    }

  fclose(in);
  fclose(out);
  return SUCCESS;
}


/* display usage */
void usage(char *command)
{
  printf("\nUsage: %s [-lndsr] <in_lndsr_cfmask> [-cmask] <out_cmask_bin_file> [-acc] <acceptable_cloud_confidence(0-3, optional)>\n\n", command);
  printf("   -lndsr  <landsat_cfmask_File> (input) in binary format\n");
  printf("   -cmask  <cloud_mask_file> (output) in binary format\n");
  printf("   -acc    <cloud confidence) ((optional, 0=none; 1=low; 2=medium; 3=high)\n\n"); 
}


int readENVIHeader(GRID_SR *sensor) 
{
  char name[MAX_STRLEN];
  char  buffer[MAX_STRLEN] = "\0";
  char  *label = NULL;
  char  *tokenptr = NULL;
  char  *seperator = "=\" \t";
  
  FILE *fp;

  sprintf(name, "%s.hdr", sensor->fileName);
  if((fp=fopen(name, "r"))==NULL) {
    printf("Open file %s error!\n", name);
    return FAILURE;
  }
  
   /* process line by line */
  while(fgets(buffer, MAX_STRLEN, fp) != NULL) {
   /* get string token */
    tokenptr = strtok(buffer, seperator);
    label=tokenptr;
 
    while(tokenptr != NULL) {
 
      tokenptr = strtok(NULL, seperator);
      if(strcmp(label, "lines") == 0)
	sscanf(tokenptr, "%d", &(sensor->nrows));
      if(strcmp(label, "samples") == 0)
	sscanf(tokenptr, "%d", &(sensor->ncols));
      /* in case label (key words) is no the first word in a line */
      label = tokenptr;
    }
  }
  fclose(fp);
  return SUCCESS;
}




/* convert decimal number to binary number */
void dec2bin(unsigned short int num, int *bitpattern, int limit)
{
 
  register int i=0;       
 
  for(i=0; i<limit; i++)
    bitpattern[i]=0;
 
  for(i=0; num>0; i++)
    {
      bitpattern[i] = (int)num & 1;
      num >>= 1;
    }
}
