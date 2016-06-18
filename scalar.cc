#include <iostream>
#include <complex>
#include <fstream>
#include <iterator>
#include <string>
#include <sstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include "geometry.hh"
#include "gamma.hh"
#include "scalar.hh"

using namespace boost::numeric::ublas;
using std::complex;
using std::size_t;

extern compressed_matrix<complex<double> > gamma5;
extern compressed_matrix<complex<double> > gamma0;
extern compressed_matrix<complex<double> > gamma1;
extern compressed_matrix<complex<double> > gamma2;
extern compressed_matrix<complex<double> > gamma3;

void scalar(double cor[],
            vector< complex<double> > &u,
            vector< complex<double> > &v,
            const long unsigned int T, const long unsigned int L) {
  
  vector< complex<double> > vtmp(12);
  for(size_t t = 0; t < T; t++) {
    cor[t] = 0.;
    for(size_t x = 0; x < L; x++) {
      for(size_t y = 0; y < L; y++) {
	for(size_t z = 0; z < L; z++) {
	  size_t ix = get_index(t, x, y, z, T, L);
          size_t mu = 0;
          size_t k = gsi(ix);
          // S^\dagger(x) \gamma5\gamma_mu times v
          for(size_t i = 0; i < 12; i+=3) {
            cor[t] += real(inner_prod(conj(project(u, range(k+i, k+i+3))), gamma_factor[5][i/3]*project(v, range(k+i*3, k+i*3+3))));
          }
	}
      }
    }

    cor[t] /= double(L*L*L);

  }
  return;
}
