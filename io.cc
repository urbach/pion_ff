/* $Id: io.c,v 1.2 2006/04/18 15:29:27 urbach Exp $ */

/****************************************************
 * IO routines:
 *
 * read_lime_gauge_field_doubleprec
 *
 * read_lime_gauge_field_singleprec
 *
 * Autor: 
 *        Carsten Urbach <urbach@ifh.de>
 *
 ****************************************************/

/*
 * Note:
 * Required version of lime: >= 1.2.3
 * n_uint64_t is a lime defined type!!
 *
 */

#define _FILE_OFFSET_BITS 64

#include"lime.h" 
#include<complex>
#include<stdio.h>
#include<string.h>
#include<iostream>
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include"dml.hh"
#include "io.hh"

using namespace boost::numeric::ublas;
using std::complex;
using std::cout;
using std::endl;

#define MAXBUF 1048576

void byte_swap(void *ptr, int nmemb);
void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_singleprec(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb);
void single2double(void * out_ptr, void * in_ptr, int nmemb);
void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb);
void double2single(void * out_ptr, void * in_ptr, int nmemb);
int big_endian();

void read_cmi(vector< complex<double> > &v, const int L, const int T, const char * filename) {
  FILE * ifs;
  float tmp[24];
  ifs = fopen(filename, "r");

  for(int x = 0; x < L; x++) {
    for(int y = 0; y < L; y++) {
      for(int z = 0; z < L; z++) {
	for(int t = 0; t < T; t++) {
	  long unsigned int ix = (t*L*L*L + x*L*L + y*L + z)*12;
	  complex<double> phase(cos(t*3.141593/T), sin(t*3.141593/T));
	  fread(tmp, 24*sizeof(float), 1, ifs);
	  for(int i = 0; i < 12; i++) {
	    v(ix+i) = phase*complex<double>(tmp[2*i], tmp[2*i+1]);
	  }
	}
      }
    }
  }
  return;
}

int read_binary_spinor_data(vector< complex<double> > &v, LimeReader * limereader, 
			    const double prec, const int ts, DML_Checksum &ans,
			    const unsigned int T, const unsigned int LX, 
			    const unsigned int LY, const unsigned int LZ) {
  int status=0;
  n_uint64_t bytes;
  double tmp[24], tmp1[24];
  DML_SiteRank rank;
  float tmp2[24];
  int words_bigendian;
  words_bigendian = big_endian();

  DML_checksum_init(&ans);
  rank = (DML_SiteRank) 0;
  
  if(prec == 32) bytes = 24*sizeof(float);
  else bytes = 24*sizeof(double);
  for(unsigned int t = 0; t < T; t++){
    if(ts > -1 && (unsigned int)abs(ts) < T) {
      t = ts;
      limeReaderSeek(limereader,(n_uint64_t) 
		     (t*LZ*LY*LX)*bytes,
		     SEEK_SET);
    }
    for(unsigned int z = 0; z < LZ; z++) {
      for(unsigned int y = 0; y < LY; y++) {
	for(unsigned int x = 0; x < LX; x++) {
	  long unsigned int ix = (t*LX*LY*LZ + x*LY*LZ + y*LZ + z)*12;
	  rank = (DML_SiteRank) (((t*LZ + z)*LY + y)*LX + x);
	  if(prec == 32) {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp2, bytes);	    
	  }
	  else {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    DML_checksum_accum(&ans,rank,(char *) tmp, bytes);
	  }
	  if(!words_bigendian) {
	    if(prec == 32) {
	      byte_swap_assign_single2double(tmp1, (float*)tmp2, 24);
	    }
	    else {
	      byte_swap_assign(tmp1, tmp, 24);
	    }
	  }
	  else {
	    if(prec == 32) {
	      single2double(tmp1, (float*)tmp2, 24);
	    }
	    else memcpy(tmp1, tmp, bytes);
	  }
	  complex<double> phase(1,0); //cos(t*3.141593/T), sin(t*3.141593/T));
	  for(unsigned int i=0; i < 12; i++) {
	    v(ix + i) = phase*complex<double>(tmp1[2*i], tmp1[2*i+1]);
	  }
	  if(status < 0 && status != LIME_EOR) {
	    return(-1);
	  }
	}
      }
    }
    if(ts > -1 && ts < (int)T) {
      t = T;
    }
  }
  printf("The final checksum is %#lx %#lx\n", (unsigned long)ans.suma, (unsigned long)ans.sumb);
  return(0);
}

int read_lime_spinor(vector< complex<double> > &v, char * filename, const int position,
		     const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ) {
  FILE * ifs;
  int status=0, getpos=-1;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  n_uint64_t prec = 32;
  DML_Checksum checksum;
  
  if((ifs = fopen(filename, "r")) == (FILE*)NULL) {
    fprintf(stderr, "Error opening file %s\n", filename);
    return(-1);
  }

  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    return(-1);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);

    if(strcmp("scidac-binary-data",header_type) == 0) getpos++;
    printf("... found record of type %s pos = %d!\n", header_type, getpos);
    if(getpos == position) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no scidac-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    return(-1);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(double))) prec = 64;
  else if(bytes == LX*LY*LZ*T*(uint64_t)(24*sizeof(float))) prec = 32;
  else {
    fprintf(stderr, "wrong length in eospinor: bytes = %lu, not %lu. Aborting read!\n", 
	    (unsigned long)bytes, (unsigned long)(LX*LY*LZ*T*(uint64_t)(24*sizeof(double))));
    return(-1);
  }
  printf("# %lu Bit precision read\n", (unsigned long)prec);

  status = read_binary_spinor_data(v, limereader, prec, ts, checksum, T, LX, LY, LZ);

  if(status < 0) {
    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
	    status, filename);
    exit(500);
  }

  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}

int read_lime_gauge_field_doubleprec(matrix< complex<double> > &config, const char * filename,
				     const int T, const int LX, const int LY, const int LZ) {
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  double tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)) {
    if(bytes != (n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)/2) {
      fprintf(stderr, "Probably wrong lattice size or precision (bytes=%lu) in file %s expected %lu\n", 
	      (unsigned long)((n_uint64_t)bytes), filename, (unsigned long)((n_uint64_t)LX*LY*LZ*T*72*(n_uint64_t)sizeof(double)));
      fprintf(stderr, "Aborting...!\n");
      fflush( stdout );
      exit(501);
    }
    else {
      fclose(ifs);
      fprintf(stderr, "single precision read!\n");
      return( read_lime_gauge_field_singleprec(config, filename, T, LX, LY, LZ) );
    }
  }

  bytes = (n_uint64_t)72*sizeof(double);

  for(unsigned t = 0; t < T; t++) {
    for(unsigned z = 0; z < LZ; z++) {
      for(unsigned y = 0; y < LY; y++) {
	for(unsigned x = 0; x < LX; x++) {
	  unsigned long int p = (((t*LX+x)*LY+y)*LZ+z)*12;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign(tmp2, tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	  }
	  unsigned long int k =0;
	  //ILDG has mu-order: x,y,z,t
	  for(unsigned int mu = 1; mu < 4; mu++) {
	    for(unsigned int i = 0; i < 3; i++) {
	      for(unsigned int j = 0; j < 3; j++) {
		config (p+mu*3+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
		k++;
	      }
	    }
	  }
 	  for(unsigned int i = 0; i < 3; i++) {
 	    for(unsigned int j = 0; j < 3; j++) {
 	      config (p+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
 	      k++; 	    
	    }
 	  }
	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int read_lime_gauge_field_singleprec(matrix< complex<double> > &config, const char * filename,
				     const int T, const int LX, const int LY, const int LZ){
  FILE * ifs;
  int t, x, y, z, status;
  n_uint64_t bytes;
  char * header_type;
  LimeReader * limereader;
  float tmp[72], tmp2[72];
  int words_bigendian;

  words_bigendian = big_endian();
  ifs = fopen(filename, "r");
  if(ifs == (FILE *)NULL) {
    fprintf(stderr, "Could not open file %s\n Aborting...\n", filename);
    exit(500);
  }
  limereader = limeCreateReader( ifs );
  if( limereader == (LimeReader *)NULL ) {
    fprintf(stderr, "Unable to open LimeReader\n");
    exit(500);
  }
  while( (status = limeReaderNextRecord(limereader)) != LIME_EOF ) {
    if(status != LIME_SUCCESS ) {
      fprintf(stderr, "limeReaderNextRecord returned error with status = %d!\n", status);
      status = LIME_EOF;
      break;
    }
    header_type = limeReaderType(limereader);
    if(!strcmp("ildg-binary-data",header_type)) break;
  }
  if(status == LIME_EOF) {
    fprintf(stderr, "no ildg-binary-data record found in file %s\n",filename);
    limeDestroyReader(limereader);
    fclose(ifs);
    exit(-2);
  }
  bytes = limeReaderBytes(limereader);
  if((int)bytes != LX*LY*LZ*T*72*sizeof(float)) {
    fprintf(stderr, "Probably wrong lattice size or precision (bytes=%d) in file %s\n", (int)bytes, filename);
    fprintf(stderr, "Aborting...!\n");
    fflush( stdout );
    exit(501);
  }

  bytes = (n_uint64_t)72*sizeof(float);
  for(t = 0; t < T; t++){
    for(z = 0; z < LZ; z++){
      for(y = 0; y < LY; y++){
	for(x = 0; x < LX; x++) {
	  int p = (((t*LX+x)*LY+y)*LZ+z)*12;
	  if(!words_bigendian) {
	    status = limeReaderReadData(tmp, &bytes, limereader);
	    byte_swap_assign_singleprec(tmp2, tmp, 72);
	  }
	  else {
	    status = limeReaderReadData(tmp2, &bytes, limereader);
	  }
	  int k =0;
	  //ILDG has mu-order: x,y,z,t
	  for(int mu = 1; mu < 4; mu++) {
	    for(int i = 0; i < 3; i++) {
	      for(int j = 0; j < 3; j++) {
		config (p+mu*3+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
		k++;
	      }
	    }
	  }
	  for(int i = 0; i < 3; i++) {
	    for(int j = 0; j < 3; j++) {
	      config (p+i, j) = complex<double> (tmp2[2*k], tmp2[2*k+1]);
	      k++;
	    }
	  }

	  if(status < 0 && status != LIME_EOR) {
	    fprintf(stderr, "LIME read error occured with status = %d while reading file %s!\n Aborting...\n", 
		    status, filename);
	    exit(500);
	  }
	}
      }
    }
  }
  limeDestroyReader(limereader);
  fclose(ifs);
  return(0);
}


int big_endian(){
  union{
    int l;
    char c[sizeof(int)];
  } u;

  u.l=1;
  return(u.c[sizeof(int) - 1] == 1);
}

void byte_swap(void * ptr, int nmemb){
  int j;
  char char_in[4];
  char * in_ptr;
  int * int_ptr;

  for(j = 0, int_ptr = (int *) ptr; j < nmemb; j++, int_ptr++) {
    in_ptr = (char *) int_ptr;
    
    char_in[0] = in_ptr[0];
    char_in[1] = in_ptr[1];
    char_in[2] = in_ptr[2];
    char_in[3] = in_ptr[3];

    in_ptr[0] = char_in[3];
    in_ptr[1] = char_in[2];
    in_ptr[2] = char_in[1];
    in_ptr[3] = char_in[0];
  }
}

void byte_swap_assign(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr, * double_out_ptr;

  double_in_ptr = (double *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) double_in_ptr;
    char_out_ptr = (char *) double_out_ptr;
    
    char_out_ptr[7] = char_in_ptr[0];
    char_out_ptr[6] = char_in_ptr[1];
    char_out_ptr[5] = char_in_ptr[2];
    char_out_ptr[4] = char_in_ptr[3];
    char_out_ptr[3] = char_in_ptr[4];
    char_out_ptr[2] = char_in_ptr[5];
    char_out_ptr[1] = char_in_ptr[6];
    char_out_ptr[0] = char_in_ptr[7];
    double_in_ptr++;
    double_out_ptr++;
  }
}

void byte_swap_assign_singleprec(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  float * float_in_ptr, * float_out_ptr;

  float_in_ptr = (float *) in_ptr;
  float_out_ptr = (float *) out_ptr;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) float_in_ptr;
    char_out_ptr = (char *) float_out_ptr;
    
    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];
    float_in_ptr++;
    float_out_ptr++;
  }
}

void single2double(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) in_ptr;
  double * double_ptr = (double*) out_ptr;

  for(i = 0; i < nmemb; i++) {
    (*double_ptr) = (double) (*float_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void double2single(void * out_ptr, void * in_ptr, int nmemb) {
  int i;
  float * float_ptr = (float*) out_ptr;
  double * double_ptr = (double*) in_ptr;

  for(i = 0; i < nmemb; i++) {
    (*float_ptr) = (float) (*double_ptr);

    float_ptr++;
    double_ptr++;
  }

}

void byte_swap_assign_single2double(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_out_ptr;
  float * float_in_ptr;
  float tmp;

  float_in_ptr = (float *) in_ptr;
  double_out_ptr = (double *) out_ptr;
  char_out_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    char_in_ptr = (char *) float_in_ptr;
    
    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];
    (*double_out_ptr) = (double) tmp;
    float_in_ptr++;
    double_out_ptr++;
  }
}

void byte_swap_assign_double2single(void * out_ptr, void * in_ptr, int nmemb){
  int j;
  char * char_in_ptr, * char_out_ptr;
  double * double_in_ptr;
  float * float_out_ptr;
  float tmp;

  float_out_ptr = (float *) out_ptr;
  double_in_ptr = (double *) in_ptr;
  char_in_ptr = (char *) &tmp;
  for(j = 0; j < nmemb; j++){
    tmp = (float) (*double_in_ptr);
    char_out_ptr = (char*) float_out_ptr;

    char_out_ptr[3] = char_in_ptr[0];
    char_out_ptr[2] = char_in_ptr[1];
    char_out_ptr[1] = char_in_ptr[2];
    char_out_ptr[0] = char_in_ptr[3];

    float_out_ptr++;
    double_in_ptr++;
  }
}

