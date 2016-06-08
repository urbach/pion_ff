#ifndef _GAMMA_HH
#define _GAMMA_HH

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/program_options.hpp>
#include "gamma.hh"

using namespace boost::numeric::ublas;
using std::complex;

void rotate_etmc_ukqcd(vector< complex<double> > &u, 
		       vector< complex<double> > &v,
		       const compressed_matrix<complex<double> > &gamma0,
		       const compressed_matrix<complex<double> > &gamma5,
		       const int vol);

void define_etmc_gammas(compressed_matrix<complex<double> > &gamma0,
			compressed_matrix<complex<double> > &gamma1,
			compressed_matrix<complex<double> > &gamma2,
			compressed_matrix<complex<double> > &gamma3,
			compressed_matrix<complex<double> > &gamma5);

void define_ukqcd_gammas(compressed_matrix<complex<double> > &gamma0,
			 compressed_matrix<complex<double> > &gamma1,
			 compressed_matrix<complex<double> > &gamma2,
			 compressed_matrix<complex<double> > &gamma3,
			 compressed_matrix<complex<double> > &gamma5);

const int gamma_permutation[6][4] = {

  {2, 3, 0, 1}, // gamma 0
  {3, 2, 1, 0}, // gamma 1
  {3, 2, 1, 0}, // gamma 2
  {2, 3, 0, 1}, // gamma 3
  {0, 1, 2, 3}, // 1
  {0, 1, 2, 3} // gamma 5
};

const complex<double> gamma_factor[6][4] = {

  {complex<double>(1.), complex<double>(1.), complex<double>(1.), complex<double>(1.)}, 
  {complex<double>(0,-1), complex<double>(0,-1), complex<double>(0,1.), complex<double>(0,1.)},
  {complex<double>(1.), complex<double>(-1.), complex<double>(-1.), complex<double>(1.)}, 
  {complex<double>(0,-1.), complex<double>(0,1.), complex<double>(0,1.), complex<double>(0,-1.)},
  {complex<double>(1.), complex<double>(1.), complex<double>(1.), complex<double>(1.)},
  {complex<double>(1.), complex<double>(1.), complex<double>(-1.), complex<double>(-1.)}
};


#endif
