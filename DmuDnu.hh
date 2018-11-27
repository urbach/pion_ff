#ifndef _DRMUDRNU_HH
#define _DRMUDRNU_HH

#include <complex>
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

void DmuDnu(double cor[],
	    vector< complex<double> > &u,
	    vector< complex<double> > &v,
	    matrix< complex<double> > &config,
	    const int T, const int L,
	    size_t kappa, size_t mu, size_t nu);

#endif
