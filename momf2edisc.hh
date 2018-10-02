#ifndef _MOMF2EDISC_HH
#define _MOMF2EDISC_HH

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>

using namespace boost::numeric::ublas;
using std::complex;

void moment_f2e_disc(double cor[],
		     vector< complex<double> > &u,
		     vector< complex<double> > &v,
		     matrix< complex<double> > &config,
		     const int T, const int L);

#endif
