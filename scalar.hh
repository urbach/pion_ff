#ifndef _SCALAR_HH
#define _SCALAR_HH

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

void scalar(double cor[],
            vector< complex<double> > &u,
            vector< complex<double> > &v,
            const long unsigned int T, const long unsigned int L);

#endif
