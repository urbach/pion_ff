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
#include "gamma.hh"

using namespace boost::numeric::ublas;
using std::complex;



void rotate_etmc_ukqcd(vector< complex<double> > &u, 
		       vector< complex<double> > &v,
		       const compressed_matrix<complex<double> > &gamma0,
		       const compressed_matrix<complex<double> > &gamma5,
		       const int vol) {
  double fac = 1./sqrt(2.);
  compressed_matrix<complex<double> > g = fac*(gamma0-gamma5);
  for(int i = 0; i < vol; i++) {
    vector_range< vector< complex<double> > > vr(v, range(12*i, 12*(i+1)));
    vector_range< vector< complex<double> > > ur(u, range(12*i, 12*(i+1)));
    ur = prod(g, vr);
  }
}

void define_etmc_gammas(compressed_matrix<complex<double> > &gamma0,
			compressed_matrix<complex<double> > &gamma1,
			compressed_matrix<complex<double> > &gamma2,
			compressed_matrix<complex<double> > &gamma3,
			compressed_matrix<complex<double> > &gamma5) {
  gamma0.clear();
  gamma1.clear();
  gamma2.clear();
  gamma3.clear();
  gamma5.clear();
  // g5 = g0 g1 g2 g3
  for(int i = 0; i < 12; i++) {
    if(i > 5) {
      gamma5(i,i) = complex<double>(-1.);
    }
    else {
      gamma5(i,i) = complex<double>(1.);
    }
  }
  
  for(int j = 0; j < 3; j++) {
    gamma0(0+j, 6+j) = complex<double>(1.);
    gamma0(3+j, 9+j) = complex<double>(1.);
    gamma0(6+j, 0+j) = complex<double>(1.);
    gamma0(9+j, 3+j) = complex<double>(1.);
  }

  for(int j = 0; j < 3; j++) {
    gamma1(0+j, 9+j) = complex<double>(0, 1.);
    gamma1(3+j, 6+j) = complex<double>(0, 1.);
    gamma1(6+j, 3+j) = complex<double>(0, -1.);
    gamma1(9+j, 0+j) = complex<double>(0, -1.);
  }
  for(int j = 0; j < 3; j++) {
    gamma2(0+j, 9+j) = complex<double>(1.);
    gamma2(3+j, 6+j) = complex<double>(-1.);
    gamma2(6+j, 3+j) = complex<double>(-1.);
    gamma2(9+j, 0+j) = complex<double>(1.);
  }
  for(int j = 0; j < 3; j++) {
    gamma3(0+j, 6+j) = complex<double>(0., 1.);
    gamma3(3+j, 9+j) = complex<double>(0., -1.);
    gamma3(6+j, 0+j) = complex<double>(0., -1.);
    gamma3(9+j, 3+j) = complex<double>(0., 1.);
  }
  return;
}

void define_ukqcd_gammas(compressed_matrix<complex<double> > &gamma0,
			compressed_matrix<complex<double> > &gamma1,
			compressed_matrix<complex<double> > &gamma2,
			compressed_matrix<complex<double> > &gamma3,
			compressed_matrix<complex<double> > &gamma5) {

  gamma0.clear();
  gamma1.clear();
  gamma2.clear();
  gamma3.clear();
  gamma5.clear();
  for(int i = 0; i < 12; i++) {
    if(i > 5) {
      gamma0(i,i) = complex<double>(-1.);
    }
    else {
      gamma0(i,i) = complex<double>(1.);
    }
  }

  // g5 = g1 g2 g3 g0  
  for(int j = 0; j < 3; j++) {
    gamma5(0+j, 6+j) = complex<double>(1.);
    gamma5(3+j, 9+j) = complex<double>(1.);
    gamma5(6+j, 0+j) = complex<double>(1.);
    gamma5(9+j, 3+j) = complex<double>(1.);
  }
  for(int j = 0; j < 3; j++) {
    gamma1(0+j, 9+j) = complex<double>(0, 1.);
    gamma1(3+j, 6+j) = complex<double>(0, 1.);
    gamma1(6+j, 3+j) = complex<double>(0, -1.);
    gamma1(9+j, 0+j) = complex<double>(0, -1.);
  }
  for(int j = 0; j < 3; j++) {
    gamma2(0+j, 9+j) = complex<double>(1.);
    gamma2(3+j, 6+j) = complex<double>(-1.);
    gamma2(6+j, 3+j) = complex<double>(-1.);
    gamma2(9+j, 0+j) = complex<double>(1.);
  }
  for(int j = 0; j < 3; j++) {
    gamma3(0+j, 6+j) = complex<double>(0., 1.);
    gamma3(3+j, 9+j) = complex<double>(0., -1.);
    gamma3(6+j, 0+j) = complex<double>(0., -1.);
    gamma3(9+j, 3+j) = complex<double>(0., 1.);
  }
  return;
}
