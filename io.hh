/* $Id: io.h,v 1.2 2006/04/18 15:29:27 urbach Exp $ */
#ifndef _IO_H
#define _IO_H

#include <complex>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"


using namespace boost::numeric::ublas;
using std::complex;

void read_cmi(vector< complex<double> > &v, const int L, const int T, const char * filename);

int read_lime_spinor(vector< complex<double> > &v, char * filename, const int position,
		     const int ts,
		     const unsigned int T, const unsigned int LX, const unsigned int LY, const unsigned int LZ);

int read_lime_gauge_field_doubleprec(matrix< std::complex<double> > &config, const char * filename,
				     const int T, const int LX, const int LY, const int LZ);

int read_lime_gauge_field_singleprec(matrix< std::complex<double> > &config, const char * filename,
				     const int T, const int LX, const int LY, const int LZ);


#endif
